#!/usr/bin/env python3
"""Fragility ensemble driver for the conductor model (standard library only).

From one base deck and a plan (wind-speed stripes x wind directions x gust seeds) this script

  1. writes one run directory per member under <output_dir>/members, each with a deck that differs
     from the base deck only in the inputs listed for it in <output_dir>/manifest.json, with
     conductors.intensity_measure = wind switched on;
  2. runs the members one after another (subcommand run), waiting for a quiet machine before each
     run and pausing the run while the machine is busy; a member that finished is not run again;
  3. gathers every member's conductors/demand_table.csv and conductors/intensity.csv into two long
     tables keyed by member, stripe, direction and seed (subcommand collect), and lists every member
     that failed, did not run or lacks a table in members_status.csv.

    python3 fragility_ensemble.py run --plan plan.json --dry-run
    python3 fragility_ensemble.py run --plan plan.json --exe <build>/Exec/erf_exec --ranks 2
    python3 fragility_ensemble.py collect --plan plan.json

How a member's wind is set is described in README.md beside this file and in
Docs/sphinx_doc/theory/Conductors.rst ("An ensemble of runs for fragility curves").
"""
import argparse
import csv
import hashlib
import json
import math
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time

TOOL = "fragility_ensemble"
MANIFEST_FORMAT = 1
MANIFEST = "manifest.json"
STATUS = "member_status.json"
RUN_LOG = "run.log"
LOCK = "run.lock"
DEMAND_ALL = "demand_table_all.csv"
INTENSITY_ALL = "intensity_all.csv"
STATUS_ALL = "members_status.csv"
KEY_COLUMNS = ("member", "stripe", "nominal_direction", "seed")
PLAN_KEYS = ("base_deck", "output_dir", "speeds", "directions", "seeds", "reference_height",
             "side_type", "inputs", "support_files")
SIDE_TYPES = ("SlipWall", "NoSlipWall", "Symmetry")
# inputs the members' wind is derived from: a plan that set them would run another wind than the manifest's
DRIVER_KEYS = ("conductors.intensity_measure", "conductors.prescribed_velocity", "erf.init_type",
               "erf.input_sounding_file", "geometry.is_periodic", "geometry.prob_lo", "geometry.prob_hi",
               "amr.n_cell", "amr.restart", "erf.restart", "erf.input_bndry_planes", "erf.abl_driver_type",
               "erf.sounding_type")
LATERAL_FACES = ("xlo", "xhi", "ylo", "yhi")
MAX_SPEED = 150.0      # m/s: above any wind a line is designed for
MAX_MEMBERS = 100000
PARALLEL = 1.0e-12     # a direction cosine below this is a wind along the face


class Stop(Exception):
    """An input the driver refuses; the message names the option."""


def stop(message):
    raise Stop(message)


def fmt(x):
    """A number as the decks and wind files carry it: 10 significant digits, no negative zero."""
    s = "%.10g" % x
    return "0" if float(s) == 0.0 else s


def read_text(path):
    """A text file as the driver reads every file: UTF-8, any other byte kept as it is."""
    with open(path, encoding="utf-8", errors="surrogateescape", newline=None) as f:
        return f.read()


def write_text(path, text):
    with open(path, "w", encoding="utf-8", errors="surrogateescape") as f:
        f.write(text)


def is_number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


# ---------------------------------------------------------------------------------------------
# The deck: AMReX ParmParse text, one "key = value" per line, FILE = <name> includes read in place
# ---------------------------------------------------------------------------------------------

def split_comment(text):
    quote = None
    for i, ch in enumerate(text):
        if quote:
            if ch == quote:
                quote = None
        elif ch == '"':
            quote = ch
        elif ch == "#":
            return text[:i], text[i:]
    return text, ""


def tokens_of(value):
    lex = shlex.shlex(value, posix=True)
    lex.whitespace_split, lex.quotes, lex.commenters, lex.escape = True, '"', "", ""
    try:
        return list(lex)
    except ValueError:
        stop("cannot read the value '%s' (unbalanced quotes)" % value)


class Deck:
    """The lines of a deck in order, each [key, value, text]; key is None for a comment or a blank."""

    def __init__(self, path):
        self.path = path
        self.directory = os.path.dirname(os.path.abspath(path))
        self.lines = []
        self._read(os.path.abspath(path), [])

    def _read(self, path, parents):
        if path in parents:
            stop("base_deck: %s includes itself through FILE" % path)
        try:
            raw = read_text(path).splitlines()
        except OSError as e:
            stop("base_deck: cannot read %s (%s)" % (path, e.strerror))
        for text in raw:
            body, _ = split_comment(text)
            if not body.strip():
                self.lines.append([None, None, text])
                continue
            key, sep, value = body.partition("=")
            key, value = key.strip(), value.strip()
            if not sep or not key or len(key.split()) != 1 or "=" in split_quoted(value):
                stop("base_deck: cannot read the line '%s' of %s: the driver needs one 'key = value' per line"
                     % (text.strip(), path))
            if key == "FILE":
                names = tokens_of(value)
                if len(names) != 1:
                    stop("base_deck: FILE needs one file name in %s" % path)
                self.lines.append([None, None, "# %s: the lines of %s" % (TOOL, names[0])])
                self._read(os.path.join(os.path.dirname(path), names[0]), parents + [path])
            else:
                self.lines.append([key, value, text])

    def get(self, key):
        """The value of the last line that sets key (the one ParmParse's query returns), or None."""
        found = None
        for k, v, _ in self.lines:
            if k == key:
                found = v
        return found

    def tokens(self, key):
        v = self.get(key)
        return None if v is None else tokens_of(v)

    def keys_under(self, prefix):
        """{sub: value} of the keys prefix.sub, in the deck's order."""
        out = {}
        for k, v, _ in self.lines:
            if k is not None and k.startswith(prefix + "."):
                out[k[len(prefix) + 1:]] = v
        return out

    def text_with(self, changes):
        """The deck's text with changes {key: value, or None to remove the key} applied: the first
        line of a changed key is rewritten in place, its other lines are dropped, and keys the deck
        does not have are added at the end."""
        out, done = [], set()
        for k, _, text in self.lines:
            if k is None or k not in changes:
                out.append(text)
            elif changes[k] is not None and k not in done:
                out.append("%s = %s" % (k, changes[k]))
                done.add(k)
        new = [k for k in changes if changes[k] is not None and k not in done]
        if new:
            out.append("")
            out.append("# set by %s" % TOOL)
            out.extend("%s = %s" % (k, changes[k]) for k in new)
        return "\n".join(out) + "\n"


def split_quoted(value):
    """value with its quoted parts removed, to look for a second '=' outside quotes."""
    out, quote = [], None
    for ch in value:
        if quote:
            if ch == quote:
                quote = None
        elif ch == '"':
            quote = ch
        else:
            out.append(ch)
    return "".join(out)


def same_value(a, b):
    """Two deck values are the same when their tokens are (quotes aside), numbers compared as numbers."""
    if a is None or b is None:
        return a is b
    ta, tb = tokens_of(a), tokens_of(b)
    if len(ta) != len(tb):
        return False
    for x, y in zip(ta, tb):
        if x == y:
            continue
        try:
            if float(x) != float(y):
                return False
        except ValueError:
            return False
    return True


# ---------------------------------------------------------------------------------------------
# The wind files: ERF's input sounding (z theta qv u v under a header line) and the inflow
# profile of <face>.dirichlet_file (z u v w [theta])
# ---------------------------------------------------------------------------------------------

class WindFile:
    def __init__(self, path, kind, option):
        self.path, self.kind = path, kind
        self.ucol = 3 if kind == "sounding" else 1
        self.rows = []      # [tokens, z, u, v] or [None, text] for a line kept as it is
        try:
            self.text = read_text(path)
            raw = self.text.splitlines()
        except OSError as e:
            stop("%s: cannot read %s (%s)" % (option, path, e.strerror))
        need = 5 if kind == "sounding" else 4
        header_left = 1 if kind == "sounding" else 0
        last_z = None
        for text in raw:
            tok = text.split()
            if not tok or tok[0].startswith("#"):
                self.rows.append([None, text])
                continue
            if header_left:
                header_left -= 1
                self.rows.append([None, text])
                continue
            try:
                z, u, v = float(tok[0]), float(tok[self.ucol]), float(tok[self.ucol + 1])
            except (ValueError, IndexError):
                stop("%s: cannot read the row '%s' of %s (%d columns expected)" % (option, text.strip(), path, need))
            if len(tok) < need:
                stop("%s: the row '%s' of %s has fewer than %d columns" % (option, text.strip(), path, need))
            if not all(math.isfinite(x) for x in (z, u, v)):
                stop("%s: the row '%s' of %s is not finite" % (option, text.strip(), path))
            if last_z is not None and z <= last_z:
                stop("%s: the heights of %s must increase" % (option, path))
            last_z = z
            self.rows.append([tok, z, u, v])
        self.data = [r for r in self.rows if r[0] is not None]
        if not self.data:
            stop("%s: %s holds no wind rows" % (option, path))

    def at(self, height):
        """(u, v) at height by linear interpolation between the rows."""
        d = self.data
        if height < d[0][1] or height > d[-1][1]:
            stop("reference_height must be in [%s, %s] m, the heights of %s; got %s"
                 % (fmt(d[0][1]), fmt(d[-1][1]), self.path, fmt(height)))
        for lo, hi in zip(d, d[1:]):
            if height <= hi[1]:
                if lo[1] < height < hi[1] and math.hypot(lo[2], lo[3]) == 0.0:
                    stop("reference_height = %s m lies between the row at %s m, where the wind is zero, and the row "
                         "at %s m of %s: the rows do not resolve the profile there, so give a height from %s m up"
                         % (fmt(height), fmt(lo[1]), fmt(hi[1]), self.path, fmt(hi[1])))
                w = (height - lo[1]) / (hi[1] - lo[1])
                return lo[2] + w * (hi[2] - lo[2]), lo[3] + w * (hi[3] - lo[3])
        return d[-1][2], d[-1][3]   # a single row at that height

    def text_with(self, factor, cosr, sinr):
        """The file with every row's horizontal wind turned by the rotation and all three
        components multiplied by factor; the other columns and lines are kept as they are."""
        if factor == 1.0 and cosr == 1.0 and sinr == 0.0:
            return self.text
        out = []
        for r in self.rows:
            if r[0] is None:
                out.append(r[1])
                continue
            tok = list(r[0])
            tok[self.ucol], tok[self.ucol + 1] = (fmt(c) for c in turned(r[2], r[3], factor, cosr, sinr))
            if self.kind == "profile":
                tok[3] = fmt(factor * float(tok[3]))
            out.append(" ".join(tok))
        return "\n".join(out) + "\n"


    def max_speed(self):
        return max(math.hypot(r[2], r[3]) for r in self.data)


def turned(u, v, factor, cosr, sinr):
    """(u, v) turned and scaled; a component that is round-off of the speed is zero."""
    a, b = factor * (cosr * u - sinr * v), factor * (sinr * u + cosr * v)
    tiny = 1.0e-12 * math.hypot(a, b)
    return (0.0 if abs(a) <= tiny else a), (0.0 if abs(b) <= tiny else b)


def rotation(degrees):
    """(cos, sin) of an angle in degrees, exact at multiples of 90."""
    quarter = degrees / 90.0
    if quarter == round(quarter):
        return ((1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0))[int(round(quarter)) % 4]
    a = math.radians(degrees)
    return math.cos(a), math.sin(a)


def direction_of(u, v):
    """Degrees from +x towards +y of the direction the wind blows towards, in (-180, 180]."""
    return math.degrees(math.atan2(v, u))


# ---------------------------------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------------------------------

def deck_value(key, x):
    """A plan value as deck text."""
    if isinstance(x, bool):
        return "true" if x else "false"
    if is_number(x):
        return repr(x)
    if isinstance(x, str) and x.strip() and not any(c in x for c in "\n#="):
        return x.strip()
    if isinstance(x, list) and x and all(is_number(i) or (isinstance(i, str) and i.strip()) for i in x):
        return " ".join(deck_value(key, i) for i in x)
    stop("inputs: the value of %s must be a number, true or false, a string (without '#' or '=') or a list of them"
         % key)


def read_inputs(where, d):
    if not isinstance(d, dict):
        stop("%s must be an object of deck inputs, {\"key\": value}" % where)
    out = {}
    for k, v in d.items():
        if not k or len(k.split()) != 1 or "=" in k or k == "FILE":
            stop("%s: '%s' is not a deck key" % (where, k))
        out[k] = deck_value(k, v)
    return out


class Plan:
    def __init__(self, path):
        self.path = os.path.abspath(path)
        self.directory = os.path.dirname(self.path)
        try:
            with open(self.path) as f:
                p = json.load(f)
        except OSError as e:
            stop("--plan: cannot read %s (%s)" % (path, e.strerror))
        except ValueError as e:
            stop("--plan: %s is not JSON (%s)" % (path, e))
        if not isinstance(p, dict):
            stop("--plan: %s must hold one JSON object" % path)
        unknown = sorted(k for k in p if k not in PLAN_KEYS)
        if unknown:
            stop("--plan: unknown key(s) %s; the plan's keys are %s" % (", ".join(unknown), ", ".join(PLAN_KEYS)))

        self.base_deck = p.get("base_deck")
        if not isinstance(self.base_deck, str) or not self.base_deck.strip():
            stop("base_deck must name the base deck (a path relative to the plan)")
        self.output_dir = p.get("output_dir", "ensemble")
        if not isinstance(self.output_dir, str) or not self.output_dir.strip():
            stop("output_dir must name a directory (a path relative to the plan)")

        speeds = p.get("speeds")
        if not isinstance(speeds, list) or not speeds:
            stop("speeds must be a list of at least one wind speed (m/s)")
        self.speeds, self.speed_inputs = [], []
        for s in speeds:
            extra = {}
            if isinstance(s, dict):
                bad = sorted(k for k in s if k not in ("speed", "inputs"))
                if bad:
                    stop("speeds: unknown key(s) %s in a stripe; a stripe is a number or {\"speed\": .., \"inputs\": {..}}"
                         % ", ".join(bad))
                extra = read_inputs("speeds: a stripe's inputs", s.get("inputs", {}))
                s = s.get("speed")
            if not is_number(s) or not 0.0 < s <= MAX_SPEED:
                stop("speeds must be in (0, %s] m/s, got %s" % (fmt(MAX_SPEED), json.dumps(s)))
            if float(s) in self.speeds:
                stop("speeds: %s m/s is listed twice" % fmt(s))
            self.speeds.append(float(s))
            self.speed_inputs.append(extra)

        self.directions = None
        if "directions" in p:
            d = p["directions"]
            if not isinstance(d, list) or not d:
                stop("directions must be a list of at least one direction (degrees)")
            self.directions = []
            for x in d:
                if not is_number(x) or not -180.0 < x <= 180.0:
                    stop("directions must be in (-180, 180] degrees (the direction the wind blows towards, from +x "
                         "towards +y), got %s" % json.dumps(x))
                if float(x) in self.directions:
                    stop("directions: %s degrees is listed twice" % fmt(x))
                self.directions.append(float(x))

        self.seeds = None
        if "seeds" in p:
            s = p["seeds"]
            if not isinstance(s, list) or not s:
                stop("seeds must be a list of at least one gust seed")
            for x in s:
                if not isinstance(x, int) or isinstance(x, bool) or not 0 <= x <= 2147483647:
                    stop("seeds must be integers in [0, 2147483647], got %s" % json.dumps(x))
            if len(set(s)) != len(s):
                stop("seeds: a seed is listed twice")
            self.seeds = list(s)

        self.reference_height = p.get("reference_height")
        if self.reference_height is not None and (not is_number(self.reference_height)
                                                  or not 0.0 < self.reference_height <= 500.0):
            stop("reference_height must be in (0, 500] m, got %s" % json.dumps(self.reference_height))

        self.side_type = p.get("side_type")
        if self.side_type is not None and self.side_type not in SIDE_TYPES:
            stop("side_type must be one of %s, got %s" % (", ".join(SIDE_TYPES), json.dumps(self.side_type)))

        self.inputs = read_inputs("inputs", p.get("inputs", {}))

        files = p.get("support_files", [])
        if not isinstance(files, list) or not all(isinstance(f, str) and f.strip() for f in files):
            stop("support_files must be a list of file names (paths relative to the base deck)")
        for f in files:
            if os.path.isabs(f) or ".." in f.replace("\\", "/").split("/"):
                stop("support_files: '%s' must lie in the base deck's directory or below it" % f)
        if len(set(files)) != len(files):
            stop("support_files: a file is listed twice")
        self.support_files = list(files)

        n = len(self.speeds) * len(self.directions or [0]) * len(self.seeds or [0])
        if n > MAX_MEMBERS:
            stop("speeds x directions x seeds gives %d members; at most %d" % (n, MAX_MEMBERS))

    def path_of(self, name):
        return os.path.normpath(os.path.join(self.directory, name))


# ---------------------------------------------------------------------------------------------
# The members
# ---------------------------------------------------------------------------------------------

def name_tag(x):
    return fmt(x).replace("-", "m").replace(".", "p")


def lateral_faces(deck):
    """The x and y faces of the dimensions that are not periodic: [(dimension, low face, high face)]."""
    per = deck.tokens("geometry.is_periodic") or ["0", "0", "0"]
    try:
        flags = [int(t) for t in per[:2]]
    except ValueError:
        stop("base_deck: cannot read geometry.is_periodic = %s" % " ".join(per))
    if len(flags) < 2:
        stop("base_deck: geometry.is_periodic needs a value per dimension")
    return [(d, lo, hi) for d, (lo, hi) in enumerate((("xlo", "xhi"), ("ylo", "yhi"))) if flags[d] == 0]


def face_type(deck, face):
    t = deck.tokens(face + ".type")
    if not t:
        stop("base_deck: %s.type is not set" % face)
    return t[0]


def digest(text):
    return hashlib.sha256(text.encode("utf-8", "surrogateescape")).hexdigest()


class Ensemble:
    """The base deck read against the plan, and the members it gives."""

    def __init__(self, plan):
        self.plan = plan
        self.deck = Deck(plan.path_of(plan.base_deck))
        deck = self.deck
        self.deck_name = os.path.basename(deck.path)
        self.diagnostics_dir = (deck.tokens("conductors.diagnostics_dir") or ["conductors"])[0]
        if os.path.isabs(self.diagnostics_dir) or ".." in self.diagnostics_dir.replace("\\", "/").split("/"):
            stop("base_deck: conductors.diagnostics_dir = %s must lie inside the run directory, or the members "
                 "would write their tables to one place" % self.diagnostics_dir)
        for where, given in [("inputs", plan.inputs)] + [("speeds: a stripe's inputs", e) for e in plan.speed_inputs]:
            for key in given:
                if key in DRIVER_KEYS + ("conductors.diagnostics_dir",):
                    stop("%s: %s cannot be given: the driver sets it or reads the base deck's" % (where, key))
        model = deck.tokens("erf.structures.model")
        if plan.inputs.get("erf.structures.model", model and model[0]) != "conductors":
            stop("base_deck: %s does not set erf.structures.model = conductors" % plan.base_deck)

        self.files = {}          # member file name -> WindFile
        self.sounding = self.profile = None
        self.faces = lateral_faces(deck)
        self.template = {}
        if deck.get("conductors.prescribed_velocity") is not None:
            self.mode = "prescribed"
            if plan.reference_height is not None:
                stop("reference_height is not used with conductors.prescribed_velocity (the wind has no profile)")
            if plan.side_type is not None:
                stop("side_type is not used with conductors.prescribed_velocity (the flow's boundaries are not changed)")
            v = deck.tokens("conductors.prescribed_velocity")
            try:
                self.prescribed = [float(t) for t in v]
            except ValueError:
                self.prescribed = []
            if len(self.prescribed) != 3:
                stop("base_deck: conductors.prescribed_velocity needs three components")
            u0, v0 = self.prescribed[0], self.prescribed[1]
            self.reference_file = None
        else:
            self.mode = "sounding"
            self._read_sounding_mode()
            height = plan.reference_height if plan.reference_height is not None else 10.0
            self.reference_height = height
            ref = self.profile if self.profile is not None else self.sounding
            self.reference_file = os.path.basename(ref.path)
            u0, v0 = ref.at(height)
        self.base_speed = math.hypot(u0, v0)
        if self.base_speed <= 0.0:
            stop("base_deck: the base wind is zero%s, so it cannot be scaled to a stripe"
                 % ("" if self.mode == "prescribed" else " at reference_height"))
        self.base_direction = direction_of(u0, v0)

        gust = plan.inputs.get("conductors.gust_type", deck.get("conductors.gust_type") or "none")
        if plan.seeds is not None and tokens_of(gust) != ["random"]:
            stop("seeds need conductors.gust_type = random in the base deck or the plan's inputs (the seed is "
                 "conductors.gust_seed); the deck has gust_type = %s" % gust)
        for name in plan.support_files:
            if not os.path.isfile(os.path.join(deck.directory, name)):
                stop("support_files: %s is not a file in the base deck's directory" % name)
        # a file the deck names that stands beside it must reach the members
        for key, value, _ in deck.lines:
            if key is None:
                continue
            for tok in tokens_of(value):
                if os.path.basename(tok) in self.files or tok in plan.support_files:
                    continue
                if not tok.replace(".", "").replace("-", "").isdigit() and os.path.isfile(os.path.join(deck.directory, tok)):
                    stop("support_files: the base deck's %s names the file %s, which the members would lack; list it "
                         "in support_files" % (key, tok))
        self.sources = {self.deck_name: digest(deck.text_with({}))}
        for name, wind in self.files.items():
            self.sources[name] = digest(wind.text)
        for name in plan.support_files:
            with open(os.path.join(deck.directory, name), "rb") as f:
                self.sources[name] = hashlib.sha256(f.read()).hexdigest()
        self.members = self._members()

    def _read_sounding_mode(self):
        deck, plan = self.deck, self.plan
        init = deck.tokens("erf.init_type") or [""]
        if init[0].lower() != "input_sounding":
            stop("base_deck: the driver sets the wind through erf.init_type = input_sounding or "
                 "conductors.prescribed_velocity; the deck has erf.init_type = %s" % (init[0] or "(not set)"))
        for key in ("erf.input_bndry_planes", "erf.abl_driver_type", "erf.sounding_type"):
            v = deck.tokens(key)
            if v and v[0].lower() not in ("0", "false", "none", "ideal"):
                stop("base_deck: %s = %s is not supported (the wind would not follow the scaled sounding)" % (key, v[0]))
        for key in ("erf.abl_pressure_grad", "erf.abl_geo_wind"):
            if deck.get(key) is not None:
                stop("base_deck: %s is not supported (a forcing the driver does not scale or turn)" % key)
        for where, given in [("inputs", plan.inputs)] + [("speeds: a stripe's inputs", e) for e in plan.speed_inputs]:
            for key in given:
                if key.split(".")[0] in LATERAL_FACES:
                    stop("%s: %s cannot be given: the driver sets the lateral faces from each member's direction"
                         % (where, key))
        for key in ("amr.restart", "erf.restart"):
            if deck.get(key) is not None:
                stop("base_deck: %s is not supported: a restart would start every member from the checkpoint's wind, "
                     "not from its stripe's" % key)
        names = deck.tokens("erf.input_sounding_file") or ["input_sounding"]
        if len(names) != 1:
            stop("base_deck: erf.input_sounding_file must name one file")
        self.sounding_key_value = names[0]
        self.sounding = WindFile(os.path.join(deck.directory, names[0]), "sounding", "base_deck")
        self.files[os.path.basename(names[0])] = self.sounding

        roles = {"inflow": [], "outflow": [], "side": []}
        for _, lo, hi in self.faces:
            for face in (lo, hi):
                t = face_type(deck, face)
                role = t.lower() if t.lower() in ("inflow", "outflow") else "side"
                if role == "side" and t.lower() not in [x.lower() for x in SIDE_TYPES]:
                    stop("base_deck: %s.type = %s is not supported; a lateral face must be Inflow, Outflow, %s or periodic"
                         % (face, t, ", ".join(SIDE_TYPES)))
                sub = deck.keys_under(face)
                sub.pop("type", None)
                roles[role].append((face, t, sub))
        for role, found in roles.items():
            for face, t, sub in found:
                for bad in ("inflow_profile", "velocity"):
                    if bad in sub:
                        stop("base_deck: %s.%s is not supported; give the inflow as %s.dirichlet_file" % (face, bad, face))
            if not found:
                continue
            first = found[0]
            for face, t, sub in found[1:]:
                if t.lower() != first[1].lower() or set(sub) != set(first[2]) \
                        or any(not same_value(sub[k], first[2][k]) for k in sub):
                    if role == "side" and plan.side_type is not None:
                        continue
                    stop("base_deck: the %s faces %s and %s carry different inputs, so the driver cannot tell which "
                         "to give a face that becomes %s%s" % (role, first[0], face, role,
                                                              "; give side_type" if role == "side" else ""))
            self.template[role] = (first[1], dict(first[2]))
        if plan.side_type is not None:
            self.template["side"] = (plan.side_type, {})
        if "inflow" in self.template:
            name = self.template["inflow"][1].get("dirichlet_file")
            if name is not None:
                name = tokens_of(name)[0]
                self.profile = WindFile(os.path.join(deck.directory, name), "profile", "base_deck")
                if os.path.basename(name) in self.files:
                    stop("base_deck: the inflow profile and the sounding have the same file name")
                self.files[os.path.basename(name)] = self.profile
                self.template["inflow"][1]["dirichlet_file"] = '"%s"' % os.path.basename(name)

    def _face_changes(self, direction, factor, cx, cy):
        """The boundary inputs of a member whose wind blows along (cx, cy): on each dimension that is
        not periodic the face the wind enters is the inflow, the one it leaves the outflow, and both
        are side faces when the wind is along them."""
        deck, changes = self.deck, {}
        for dim, lo, hi in self.faces:
            c = (cx, cy)[dim]
            if abs(c) < PARALLEL:
                want = {lo: "side", hi: "side"}
            elif c > 0.0:
                want = {lo: "inflow", hi: "outflow"}
            else:
                want = {lo: "outflow", hi: "inflow"}
            for face in (lo, hi):
                role = want[face]
                now_type = face_type(deck, face)
                now_role = now_type.lower() if now_type.lower() in ("inflow", "outflow") else "side"
                if role not in self.template:
                    if role == "side":
                        stop("side_type: the direction %s degrees blows along the %s and %s faces, and the base deck "
                             "has no side face to copy; give side_type" % (fmt(direction), lo, hi))
                    stop("directions: the direction %s degrees blows through the %s and %s faces, which the base deck "
                         "closes (%s.type = %s); the base deck needs an Inflow and an Outflow face to copy"
                         % (fmt(direction), lo, hi, face, now_type))
                new_type, new_sub = self.template[role]
                new_sub = dict(new_sub)
                if role == "inflow" and "KE" in new_sub:
                    try:
                        new_sub["KE"] = fmt(float(tokens_of(new_sub["KE"])[0]) * factor * factor)
                    except (ValueError, IndexError):
                        stop("base_deck: cannot read the inflow face's KE = %s" % new_sub["KE"])
                if role != now_role or (role == "side" and new_type.lower() != now_type.lower()):
                    changes[face + ".type"] = '"%s"' % new_type
                old_sub = deck.keys_under(face)
                old_sub.pop("type", None)
                for sub in old_sub:
                    if sub not in new_sub:
                        changes["%s.%s" % (face, sub)] = None
                for sub, value in new_sub.items():
                    if not same_value(old_sub.get(sub), value):
                        changes["%s.%s" % (face, sub)] = value
        return changes

    def _members(self):
        plan, deck = self.plan, self.deck
        directions = plan.directions if plan.directions is not None else [self.base_direction]
        seeds = plan.seeds if plan.seeds is not None else [None]
        total = len(plan.speeds) * len(directions) * len(seeds)
        width = max(3, len(str(total)))
        members = []
        for speed, extra in zip(plan.speeds, plan.speed_inputs):
            for direction in directions:
                for seed in seeds:
                    factor = speed / self.base_speed
                    turn = direction - self.base_direction if plan.directions is not None else 0.0
                    cosr, sinr = rotation(turn)
                    cx, cy = rotation(direction) if plan.directions is not None else \
                        (math.cos(math.radians(direction)), math.sin(math.radians(direction)))
                    own = {}
                    if self.mode == "prescribed":
                        u, v, w = self.prescribed
                        own["conductors.prescribed_velocity"] = "%s %s %s" % (
                            tuple(fmt(c) for c in turned(u, v, factor, cosr, sinr)) + (fmt(factor * w),))
                    else:
                        base = os.path.basename(self.sounding_key_value)
                        if base != self.sounding_key_value:
                            own["erf.input_sounding_file"] = '"%s"' % base
                        own.update(self._face_changes(direction, factor, cx, cy))
                    top = None
                    if self.mode == "sounding":
                        top = factor * max(f.max_speed() for f in self.files.values())
                        if top > MAX_SPEED:
                            stop("speeds: the stripe %s m/s at %s m makes the profile's largest wind %s m/s; the wind "
                                 "must stay within (0, %s] m/s at every height" % (
                                     fmt(speed), fmt(self.reference_height), fmt(top), fmt(MAX_SPEED)))
                    own["conductors.intensity_measure"] = "wind"
                    if seed is not None:
                        own["conductors.gust_seed"] = str(seed)
                    for where, given in (("inputs", plan.inputs), ("speeds: a stripe's inputs", extra)):
                        for key in given:
                            if key in own:
                                stop("%s: %s is set by the driver for every member and cannot be given" % (where, key))
                    changes = dict(own)
                    changes.update(plan.inputs)
                    changes.update(extra)
                    changes = {k: v for k, v in changes.items() if not same_value(deck.get(k), v)}
                    index = len(members) + 1
                    name = "m%0*d_s%s_d%s" % (width, index, name_tag(speed), name_tag(direction))
                    if seed is not None:
                        name += "_seed%d" % seed
                    members.append({
                        "member": name, "stripe": speed, "nominal_direction": direction, "seed": seed,
                        "directory": "members/" + name, "scale": factor, "rotation": turn,
                        "max_wind": top, "courant": self._courant(factor, cosr, sinr, changes),
                        "inputs": changes,
                        "_wind": (factor, cosr, sinr)})
        return members

    def _courant(self, factor, cosr, sinr, changes):
        """The largest (|u|/dx + |v|/dy) dt of the member's profile on level 0 with erf.fixed_dt, or None when
        the deck does not give them (a prescribed wind does not move the flow)."""
        if self.mode != "sounding":
            return None
        try:
            lo = [float(t) for t in self.deck.tokens("geometry.prob_lo")]
            hi = [float(t) for t in self.deck.tokens("geometry.prob_hi")]
            n = [int(t) for t in self.deck.tokens("amr.n_cell")]
            dt = float(tokens_of(changes.get("erf.fixed_dt") or self.deck.get("erf.fixed_dt"))[0])
            dx, dy = (hi[0] - lo[0]) / n[0], (hi[1] - lo[1]) / n[1]
            worst = 0.0
            for f in self.files.values():
                for r in f.data:
                    u, v = turned(r[2], r[3], factor, cosr, sinr)
                    worst = max(worst, (abs(u) / dx + abs(v) / dy) * dt)
            return float(fmt(worst)) if dt > 0.0 else None
        except (TypeError, ValueError, IndexError, ZeroDivisionError, AttributeError):
            return None

    def manifest(self):
        wind = {"mode": self.mode, "base_speed": self.base_speed, "base_direction": self.base_direction}
        if self.mode == "sounding":
            wind.update({"reference_height": self.reference_height, "reference_file": self.reference_file,
                         "files": sorted(self.files)})
        return {
            "tool": TOOL, "format": MANIFEST_FORMAT,
            "base_deck": self.plan.base_deck, "deck": self.deck_name,
            "diagnostics_dir": self.diagnostics_dir,
            "wind": wind, "support_files": self.plan.support_files, "sources": self.sources,
            "members": [{k: v for k, v in m.items() if not k.startswith("_")} for m in self.members]}

    def write(self, out):
        """Write the manifest and every member directory that is not there yet. An existing manifest
        must be the one this plan gives: a changed plan would mix two ensembles in one directory."""
        text = json.dumps(self.manifest(), indent=1) + "\n"
        path = os.path.join(out, MANIFEST)
        if os.path.exists(path):
            if read_text(path) != text:
                stop("output_dir: %s was written from another plan, base deck, wind file or support file; use a new "
                     "output_dir, or remove that directory" % path)
        elif os.path.exists(os.path.join(out, "members")):
            stop("output_dir: %s holds member directories but no %s, so they cannot be checked against the plan; use "
                 "a new output_dir, or remove that directory" % (out, MANIFEST))
        else:
            write_text(path + ".tmp", text)
            os.replace(path + ".tmp", path)     # the manifest first: every directory written is one it lists
        written = 0
        for m in self.members:
            d = os.path.join(out, m["directory"])
            deck_path = os.path.join(d, self.deck_name)
            if os.path.exists(deck_path):
                continue
            os.makedirs(d, exist_ok=True)
            try:
                for name in self.plan.support_files:
                    target = os.path.join(d, name)
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    shutil.copyfile(os.path.join(self.deck.directory, name), target)
                for name, wind in self.files.items():
                    write_text(os.path.join(d, name), wind.text_with(*m["_wind"]))
                write_text(deck_path + ".tmp", self.deck.text_with(m["inputs"]))
                os.replace(deck_path + ".tmp", deck_path)    # the deck last: its presence marks a whole directory
            except OSError as e:
                stop("output_dir: cannot write %s (%s)" % (d, e.strerror))
            written += 1
        return written

    def notes(self):
        """What a user should know before the runs start."""
        out = []
        if self.mode == "sounding" and "inflow" in self.template:
            sub = self.template["inflow"][1]
            if "KE" in sub and "density" not in sub:
                out.append("NOTE the inflow face's KE is scaled with the stripe, but ERF uses it only together with "
                           "<face>.density; without it the inflow k is zero-gradient")
        for m in self.members:
            if m["courant"] is not None and m["courant"] > 1.0:
                out.append("NOTE %s: (|u|/dx + |v|/dy) erf.fixed_dt reaches %s on level 0 (largest wind %s m/s); give "
                           "the stripe its own erf.fixed_dt" % (m["member"], fmt(m["courant"]), fmt(m["max_wind"])))
        return out


# ---------------------------------------------------------------------------------------------
# Running the members, one at a time, behind the load gate
# ---------------------------------------------------------------------------------------------

def descendants(pid):
    """The processes started by pid, children first; None when ps cannot be run."""
    try:
        out = subprocess.run(["ps", "-axo", "pid=,ppid="], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    children = {}
    for line in out.splitlines():
        p = line.split()
        if len(p) == 2 and p[0].isdigit() and p[1].isdigit():
            children.setdefault(int(p[1]), []).append(int(p[0]))
    found, queue = [], [pid]
    while queue:
        for c in children.get(queue.pop(0), []):
            found.append(c)
            queue.append(c)
    return found


def signal_all(pids, sig):
    for p in pids:
        try:
            os.kill(p, sig)
        except OSError:
            pass


def machine_load():
    try:
        return os.getloadavg()[0]
    except (OSError, AttributeError):
        stop("--no-load-gate is needed: this system does not give its load average")


def end_run(proc):
    """End a run that is not to go on: its launcher is asked to stop (SIGTERM), then every process it
    started is killed. Open MPI puts each rank in its own process group, so the processes are found
    by their parents, before the launcher goes."""
    family = descendants(proc.pid) or []
    signal_all([proc.pid] + family, signal.SIGCONT)
    proc.terminate()
    try:
        proc.wait(timeout=10)
    except subprocess.TimeoutExpired:
        pass
    signal_all(family + (descendants(proc.pid) or []), signal.SIGKILL)
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except OSError:
        pass
    proc.wait()


class Gate:
    """Wait for a quiet machine before a run and pause the run while the machine is busy."""

    def __init__(self, start, pause, resume, poll, enabled=True, load=None, log=print):
        self.start, self.pause, self.resume, self.poll, self.enabled = start, pause, resume, poll, enabled
        self.load = load or machine_load
        self.log = log

    def wait_to_start(self, name):
        if not self.enabled:
            return
        told = False
        while self.load() >= self.start:
            if not told:
                self.log("%s: waiting for the 1-minute load to fall below %s" % (name, fmt(self.start)))
                told = True
            time.sleep(self.poll)

    def watch(self, name, proc):
        """Until proc ends: stop its processes when the load passes pause, continue them below resume."""
        stopped, paused = [], False
        try:
            while True:
                try:
                    proc.wait(timeout=self.poll)
                    break
                except subprocess.TimeoutExpired:
                    pass
                if not self.enabled:
                    continue
                load = self.load()
                if not paused and load > self.pause:
                    family = descendants(proc.pid)
                    if family is None:
                        self.log("%s: load %s above %s, but the run cannot be paused: ps did not list its processes"
                                 % (name, fmt(load), fmt(self.pause)))
                        self.enabled = False
                        continue
                    stopped = family or [proc.pid]
                    signal_all(stopped, signal.SIGSTOP)
                    paused = True
                    self.log("%s: paused, load %s above %s" % (name, fmt(load), fmt(self.pause)))
                elif paused and load < self.resume:
                    signal_all(reversed(stopped), signal.SIGCONT)
                    stopped, paused = [], False
                    self.log("%s: resumed, load %s below %s" % (name, fmt(load), fmt(self.resume)))
                elif paused:
                    late = [p for p in descendants(proc.pid) or [] if p not in stopped]   # started after the pause
                    signal_all(late, signal.SIGSTOP)
                    stopped += late
        finally:
            signal_all(reversed(stopped), signal.SIGCONT)   # never leave a run stopped behind


def read_status(directory):
    try:
        s = json.loads(read_text(os.path.join(directory, STATUS)))
        return s if isinstance(s, dict) else None
    except (OSError, ValueError):
        return None


def run_members(out, manifest, command, ranks, gate, log=print):
    """Run every member that has not finished. Returns the names of the members that failed."""
    failed = []
    diag = manifest["diagnostics_dir"]
    lock = os.path.join(out, LOCK)
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        stop("output_dir: %s exists: another driver is running these members (its process id is in the file); remove "
             "the file if that driver is gone" % lock)
    os.write(fd, ("%d\n" % os.getpid()).encode())
    os.close(fd)

    def ended(signum, frame):
        raise SystemExit(128 + signum)

    old = {sig: signal.signal(sig, ended) for sig in (signal.SIGTERM, signal.SIGHUP)}
    try:
        return _run_members(out, manifest, command, ranks, gate, log, failed, diag)
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
        os.remove(lock)


def _run_members(out, manifest, command, ranks, gate, log, failed, diag):
    for m in manifest["members"]:
        name, d = m["member"], os.path.join(out, m["directory"])
        status = read_status(d)
        if status is not None and status.get("status") == "done":
            log("%s: finished earlier, skipped" % name)
            continue
        # a run starts from a clean record: an earlier attempt's status and tables must not be taken for this one's
        for stale in (STATUS, os.path.join(diag, "demand_table.csv"), os.path.join(diag, "intensity.csv")):
            if os.path.exists(os.path.join(d, stale)):
                os.remove(os.path.join(d, stale))
        gate.wait_to_start(name)
        log("%s: started (stripe %s m/s, direction %s degrees%s)" % (
            name, fmt(m["stripe"]), fmt(m["nominal_direction"]),
            "" if m["seed"] is None else ", seed %d" % m["seed"]))
        t0 = time.time()
        with open(os.path.join(d, RUN_LOG), "w") as logf:
            try:
                proc = subprocess.Popen(command + [manifest["deck"]], cwd=d, stdin=subprocess.DEVNULL, stdout=logf,
                                        stderr=subprocess.STDOUT, start_new_session=True)
            except OSError as e:
                stop("--exe/--mpiexec: cannot start '%s' (%s)" % (" ".join(command), e.strerror))
            try:
                gate.watch(name, proc)
            except BaseException:
                end_run(proc)   # interrupted: the member has no status and runs again next time
                raise
        code = proc.returncode
        state = "done" if code == 0 else "failed"
        with open(os.path.join(d, STATUS + ".tmp"), "w") as f:
            json.dump({"status": state, "returncode": code, "ranks": ranks,
                       "wall_seconds": round(time.time() - t0, 1)}, f)
            f.write("\n")
        os.replace(os.path.join(d, STATUS + ".tmp"), os.path.join(d, STATUS))
        log("%s: %s (return code %d)" % (name, state, code))
        if code != 0:
            failed.append(name)
    return failed


# ---------------------------------------------------------------------------------------------
# Collecting the members' tables
# ---------------------------------------------------------------------------------------------

def read_table(path, need):
    """(header, rows, problem) of a CSV table; problem is None for a table that can be used."""
    try:
        with open(path, newline="") as f:
            rows = [r for r in csv.reader(f) if r]
    except OSError:
        return None, None, "missing"
    if not rows:
        return None, None, "empty"
    header, body = rows[0], rows[1:]
    lacking = [c for c in need if c not in header]
    if lacking:
        return None, None, "no column %s" % ", ".join(lacking)
    clash = [c for c in KEY_COLUMNS if c in header]
    if clash:
        return None, None, "has the key column %s" % ", ".join(clash)
    if not body:
        return None, None, "no rows"
    if any(len(r) != len(header) for r in body):
        return None, None, "a row's length differs from the header's"
    return header, body, None


def log_detail(path):
    """The abort message of a run's log, or its last line."""
    try:
        with open(path, errors="replace") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
    except OSError:
        return ""
    for ln in lines:
        if "amrex::Abort" in ln or "Assertion" in ln:
            return ln[:300]
    return lines[-1][:300] if lines else ""


def collect(out, log=print):
    """Write the two long tables and the members' status. Returns the number of members not collected."""
    try:
        with open(os.path.join(out, MANIFEST)) as f:
            manifest = json.load(f)
        members, diag = manifest["members"], manifest["diagnostics_dir"]
        if manifest.get("tool") != TOOL or manifest.get("format") != MANIFEST_FORMAT:
            raise KeyError("format")
    except (OSError, ValueError, KeyError, TypeError):
        stop("output_dir: %s holds no %s written by this driver" % (out, MANIFEST))
    tables = {"demand": [None, []], "intensity": [None, []]}
    first = {"demand": None, "intensity": None}
    status_rows = []
    for m in members:
        d = os.path.join(out, m["directory"])
        key = [m["member"], fmt(m["stripe"]), fmt(m["nominal_direction"]), "" if m["seed"] is None else str(m["seed"])]
        status = read_status(d)
        code = "" if status is None or status.get("returncode") is None else str(status["returncode"])
        counts = ["", ""]
        if status is None:
            state, detail = "not_run", ""
        elif status.get("status") != "done":
            state, detail = "failed", log_detail(os.path.join(d, RUN_LOG))
        else:
            state, detail, got = "complete", "", {}
            for which, name, need in (("demand", "demand_table.csv", ("kind", "component")),
                                      ("intensity", "intensity.csv", ("point",))):
                header, body, problem = read_table(os.path.join(d, diag, name), need)
                if problem is None and tables[which][0] is not None and header != tables[which][0]:
                    problem = "its header differs from %s's" % first[which]
                if problem is not None:
                    state = ("no_%s_table" if problem == "missing" else "bad_%s_table") % which
                    detail = "%s/%s: %s" % (diag, name, problem)
                    break
                got[which] = (header, body)
            if state == "complete":     # a member goes into both tables or into neither
                for which, (header, body) in got.items():
                    if tables[which][0] is None:
                        tables[which][0], first[which] = header, m["member"]
                    tables[which][1].extend(key + r for r in body)
                counts = [str(len(got["demand"][1])), str(len(got["intensity"][1]))]
        status_rows.append(key + [state, code] + counts + [detail])

    for which, name in (("demand", DEMAND_ALL), ("intensity", INTENSITY_ALL)):
        header, rows = tables[which]
        with open(os.path.join(out, name), "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(list(KEY_COLUMNS) + (header or []))
            w.writerows(rows)
    with open(os.path.join(out, STATUS_ALL), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(list(KEY_COLUMNS) + ["status", "returncode", "demand_rows", "intensity_rows", "detail"])
        w.writerows(status_rows)
    missing = [r for r in status_rows if r[4] != "complete"]
    log("collected %d of %d members into %s and %s" % (len(members) - len(missing), len(members),
                                                       DEMAND_ALL, INTENSITY_ALL))
    for r in missing:
        log("NOT COLLECTED %s: %s%s" % (r[0], r[4], (" (" + r[8] + ")") if r[8] else ""))
    return len(missing)


# ---------------------------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------------------------

def parse(argv):
    ap = argparse.ArgumentParser(prog=TOOL + ".py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command")
    run = sub.add_parser("run", help="write the members' decks and the manifest, run the members, collect their tables")
    run.add_argument("--plan", required=True, help="the plan (JSON)")
    run.add_argument("--dry-run", action="store_true", help="write the decks and the manifest and run nothing")
    run.add_argument("--exe", help="the ERF executable (needed unless --dry-run)")
    run.add_argument("--mpiexec", default="mpiexec", help="the MPI launcher; '' runs the executable directly "
                                                          "(default mpiexec)")
    run.add_argument("--mpiexec-args", default="", help="further launcher arguments in one string, "
                     "written with '=': --mpiexec-args='--oversubscribe'")
    run.add_argument("--ranks", type=int, default=4, help="MPI ranks of each run (default 4)")
    run.add_argument("--max-ranks", type=int, default=4, help="the most ranks --ranks may ask for (default 4)")
    run.add_argument("--load-start", type=float, default=6.0,
                     help="a run starts when the 1-minute load is below this (default 6)")
    run.add_argument("--load-pause", type=float, default=12.0,
                     help="a run is paused when the load passes this (default 12)")
    run.add_argument("--load-resume", type=float, default=7.0,
                     help="a paused run goes on when the load falls below this (default 7)")
    run.add_argument("--load-poll", type=float, default=20.0, help="seconds between looks at the load (default 20)")
    run.add_argument("--no-load-gate", action="store_true", help="run without looking at the load")
    col = sub.add_parser("collect", help="gather the members' tables again")
    col.add_argument("--plan", required=True, help="the plan (JSON)")
    a = ap.parse_args(argv)
    if a.command is None:
        ap.error("give a subcommand: run or collect")
    return a


def check_run_options(a):
    """The range of every run option; returns the command that starts a member, without the deck."""
    if a.max_ranks < 1 or a.max_ranks > 1048576:
        stop("--max-ranks must be in [1, 1048576], got %d" % a.max_ranks)
    if not 1 <= a.ranks <= a.max_ranks:
        stop("--ranks must be in [1, --max-ranks = %d], got %d" % (a.max_ranks, a.ranks))
    for name, value in (("--load-start", a.load_start), ("--load-pause", a.load_pause), ("--load-resume", a.load_resume)):
        if not (math.isfinite(value) and 0.0 < value <= 100000.0):
            stop("%s must be in (0, 100000], got %s" % (name, value))
    if not a.load_resume < a.load_pause:
        stop("--load-resume must be below --load-pause = %s, got %s" % (fmt(a.load_pause), fmt(a.load_resume)))
    if not a.load_start <= a.load_pause:
        stop("--load-start must not be above --load-pause = %s, got %s" % (fmt(a.load_pause), fmt(a.load_start)))
    if not (math.isfinite(a.load_poll) and 0.0 < a.load_poll <= 3600.0):
        stop("--load-poll must be in (0, 3600] s, got %s" % a.load_poll)
    if not a.exe:
        if a.dry_run:
            return None
        stop("--exe is needed to run the members (or give --dry-run)")
    exe = os.path.abspath(a.exe)
    if not (os.path.isfile(exe) and os.access(exe, os.X_OK)):
        stop("--exe: %s is not an executable file" % a.exe)
    try:
        extra = shlex.split(a.mpiexec_args)
    except ValueError:
        stop("--mpiexec-args: cannot read '%s' (unbalanced quotes)" % a.mpiexec_args)
    if not a.mpiexec:
        if a.ranks != 1:
            stop("--ranks must be 1 without a launcher (--mpiexec ''), got %d" % a.ranks)
        if extra:
            stop("--mpiexec-args needs a launcher (--mpiexec)")
        return [exe]
    if shutil.which(a.mpiexec) is None:
        stop("--mpiexec: '%s' is not on the PATH" % a.mpiexec)
    return [a.mpiexec, "-n", str(a.ranks)] + extra + [exe]


def main(argv=None, gate=None):
    try:
        a = parse(sys.argv[1:] if argv is None else argv)
        plan = Plan(a.plan)
        out = plan.path_of(plan.output_dir)
        if a.command == "collect":
            return 1 if collect(out) else 0
        command = check_run_options(a)
        ensemble = Ensemble(plan)
        os.makedirs(out, exist_ok=True)
        written = ensemble.write(out)
        print("%d members in %s (%d directories written now); the base wind is %s m/s towards %s degrees%s" % (
            len(ensemble.members), out, written, fmt(ensemble.base_speed), fmt(ensemble.base_direction),
            "" if ensemble.mode == "prescribed" else " at %s m in %s" % (fmt(ensemble.reference_height),
                                                                        ensemble.reference_file)))
        for note in ensemble.notes():
            print(note)
        if a.dry_run:
            return 0
        if gate is None:
            gate = Gate(a.load_start, a.load_pause, a.load_resume, a.load_poll, enabled=not a.no_load_gate)
        failed = run_members(out, ensemble.manifest(), command, a.ranks, gate)
        missing = collect(out)
        return 1 if failed or missing else 0
    except Stop as e:
        print("%s: %s" % (TOOL, e), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
