#!/usr/bin/env python3
"""Validate a pheniqs-style schema JSON for the demuxolotl decoder.

Mirrors the checks in demux/decoder.go (loadConfig) and also flags inputs that
would make decode() panic or silently misbehave. Stdlib only.

Exit status: 0 = no errors (warnings allowed), 1 = errors, 2 = unreadable file.
"""
import json
import math
import re
from pathlib import Path
from harpy.common.printing import HarpyPrint

MAX_SEGMENT = 64  # decoder.sc is [64][4]float64
BASES = set("ACGT")  # anything else maps to code 4 and indexes out of range in decode()
TOKEN_RE = re.compile(r"^0:(\d+):(\d*)$")  # Go only accepts read 0, "start:end", end optional
hp = HarpyPrint()

def _no_dup_keys(pairs):
    """Reject duplicate JSON keys; Go silently keeps the last one."""
    d = {}
    for k, v in pairs:
        if k in d:
            raise ValueError(f"duplicate JSON key {k!r}")
        d[k] = v
    return d


def get(obj, key):
    """Case-insensitive lookup"""
    if not isinstance(obj, dict):
        return None
    for k, v in obj.items():
        if k.lower() == key.lower():
            return v
    return None


def is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def token_list(transform, where, err):
    """Return the list of token strings from transform.token, or None."""
    tok = get(transform, "token")
    if not isinstance(tok, list) or not all(isinstance(t, str) for t in tok):
        err(f"{where}.transform.token must be a list of strings")
        return None
    return tok


def parse_token(t, where, err):
    m = TOKEN_RE.match(t)
    if not m:
        err(f'{where}: unsupported token {t!r}, expected "0:start:end"')
        return None
    return int(m.group(1)), (int(m.group(2)) if m.group(2) else 0)


def validate(cfg):
    errors, warnings = [], []
    err, warn = errors.append, warnings.append

    if not isinstance(cfg, dict):
        return ["top level must be a JSON object"], []

    # template: decoder reads the first token's start as where the read body begins
    tpl = get(get(cfg, "template"), "transform")
    toks = token_list(tpl, "template", err) if tpl is not None else None
    if tpl is None:
        err("config has no template transform token")
    elif toks is not None:
        if not toks:
            err("config has no template transform token")
        else:
            parse_token(toks[0], "template token[0]", err)

    # molecular: exactly 4 decoders, in D, A, B, C order
    mol = get(cfg, "molecular")
    if not isinstance(mol, list):
        err('"molecular" must be a list')
        return errors, warnings
    if len(mol) != 4:
        err(f"expected 4 molecular decoders (D,A,B,C), got {len(mol)}")

    decoders = get(cfg, "decoder")
    if not isinstance(decoders, dict):
        err('"decoder" must be an object')
        decoders = {}

    segments = []  # (name, start, end) for overlap check
    for i, m in enumerate(mol):
        w = f"molecular[{i}]"
        if not isinstance(m, dict):
            err(f"{w} must be an object")
            continue

        # noise: fraction of reads expected to be junk; 1 would zero every prior
        noise = get(m, "noise")
        if not is_num(noise) or not 0 <= noise < 1:
            err(f"{w}.noise must be a number in [0, 1), got {noise!r}")
        # confidence: posterior cutoff; outside [0, 1] passes or fails every read
        conf = get(m, "confidence_threshold")
        if not is_num(conf) or not 0 <= conf <= 1:
            err(f"{w}.confidence_threshold must be a number in [0, 1], got {conf!r}")

        # token: exactly one, with a real start<end window that fits the decoder
        mt = token_list(get(m, "transform"), w, err)
        start = end = None
        if mt is not None:
            if len(mt) != 1:
                err(f"{w}: each molecular decoder needs exactly one token")
            else:
                se = parse_token(mt[0], w, err)
                if se:
                    start, end = se
                    if end <= start:
                        err(f"{w}: segment {mt[0]!r} must have end > start (end is required)")
                        start = end = None
                    elif end - start > MAX_SEGMENT:
                        err(f"{w}: barcode segment longer than {MAX_SEGMENT} bases")
                        start = end = None
                    else:
                        segments.append((w, start, end))

        # base: names the entry in "decoder" holding the barcode list
        base = get(m, "base")
        if not isinstance(base, str):
            err(f'{w}.base must be a string naming an entry in "decoder"')
            continue
        codec = get(get(decoders, base), "codec") if base in decoders else None
        if codec is None:
            err(f'{w}: decoder base {base!r} has no "decoder.{base}.codec"')
            continue
        if not isinstance(codec, dict):
            err(f"decoder.{base}.codec must be an object")
            continue

        # barcodes: used = first entry of each codec; must be ACGT, right length, unique
        seqs, ids = {}, {}
        for k, c in codec.items():
            bcs = get(c, "barcode")
            if not isinstance(bcs, list) or not all(isinstance(b, str) for b in bcs):
                err(f'decoder.{base}.codec[{k!r}].barcode must be a list of strings')
                continue
            if not bcs:
                continue  # decoder skips codecs with no barcode
            if len(bcs) > 1:
                warn(f"decoder.{base}.codec[{k!r}]: only barcode[0] is used, {len(bcs) - 1} ignored")
            bc = bcs[0].upper()
            if end is not None and len(bc) != end - start:
                err(f"barcode {k} ({bc}) is not {end - start} bases")
            bad = set(bc) - BASES
            if bad:
                err(f"barcode {k} ({bc}) has non-ACGT characters {sorted(bad)}; decoder would panic")
            if bc in seqs:
                err(f"barcode {k} duplicates {seqs[bc]} ({bc}); id lookup would be ambiguous")
            seqs[bc] = k
            short = k.removeprefix("@")
            if short in ids:
                err(f"codec ids {ids[short]!r} and {k!r} collide once '@' is stripped")
            ids[short] = k
        if not seqs:
            err(f"no barcodes for decoder base {base!r}")
            continue

        # near-identical barcodes make confident calls impossible
        uniq = [s for s in seqs if len(s) == len(next(iter(seqs)))]
        for a in range(len(uniq)):
            for b in range(a + 1, len(uniq)):
                d = sum(x != y for x, y in zip(uniq[a], uniq[b]))
                if d <= 1:
                    warn(f"decoder.{base}: {seqs[uniq[a]]} and {seqs[uniq[b]]} differ by only {d} base(s)")

    # segments in a padded read must not collide
    segments.sort(key=lambda s: s[1])
    for (n1, s1, e1), (n2, s2, e2) in zip(segments, segments[1:]):
        if s2 < e1:
            warn(f"{n1} [{s1}:{e1}) overlaps {n2} [{s2}:{e2})")
    return errors, warnings

def pheniqsSchema(path: str):
    path = Path(path)
    try:
        cfg = json.loads(path.read_text(), object_pairs_hook=_no_dup_keys)
    except (OSError, ValueError) as e:
        hp.error(
            'Failed to read decoder schema',
            f"{path.name}: {e}"
        )
    errors, warnings = validate(cfg)
    if warnings:
        hp.notice(
            "Several warnings from validating the json schema:\n- " + "\n- ".join(warnings)
        )

    if errors:
        hp.error(
            "Deconder schema failed validation",
            f"There were {len(errors)} errors detected in [yellow]{path.name}[/], the decoder will fail to run or give unreliable output until fixed.",
            "Please fix the error(s) describes below and try again.",
            "Errors in schema file",
            '- ' + "\n- ".join(errors)
        )

