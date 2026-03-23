import gzip
import itertools
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class RepeatMaskerLocus:
    """Collapsed RepeatMasker annotation on a single query sequence."""

    query_name: str
    query_start: int
    query_end: int
    family: str
    repeat_class: str
    strand: str
    identity: float
    covered_bp: int
    consensus_start: int
    consensus_end: int
    rm_id: int | None = None
    fragments: list["RepeatMaskerLocus"] = field(default_factory=list, repr=False)

    @property
    def length(self) -> int:
        return self.query_end - self.query_start

    @property
    def is_l1(self) -> bool:
        return self.family.upper().startswith("L1") or self.repeat_class.startswith("LINE/L1")


def _is_repeatmasker_data(parts: list[str]) -> bool:
    return (
        len(parts) >= 14
        and parts[0].replace(".", "", 1).isdigit()
        and parts[5].isdigit()
        and parts[6].isdigit()
    )


def _parse_consensus_interval(parts: list[str]) -> tuple[int, int]:
    nums = []
    for token in parts[11:14]:
        if token.startswith("("):
            continue
        raw = token.strip("()")
        if raw.isdigit():
            nums.append(int(raw))
    if len(nums) < 2:
        raise ValueError("Malformed RepeatMasker consensus interval")
    start = min(nums) - 1
    end = max(nums)
    if end <= start:
        raise ValueError("Malformed RepeatMasker consensus interval")
    return start, end


def _parse_rm_id(parts: list[str]) -> int | None:
    for token in reversed(parts[14:]):
        if token.isdigit():
            return int(token)
    return None


def read_repeatmasker_out(
    input_path: str | Path,
    *,
    log_path: str | Path | None = None,
) -> list[RepeatMaskerLocus]:
    """Return RepeatMasker ``.out`` fragments parsed into query-space loci."""

    opener = gzip.open if str(input_path).endswith(".gz") else open
    fragments: list[RepeatMaskerLocus] = []
    skipped: list[str] = []

    with opener(input_path, "rt") as fin:
        first_lines = [fin.readline() for _ in range(4)]
        if any("SW" in line and "perc" in line for line in first_lines):
            iterator = fin
        else:
            iterator = itertools.chain(first_lines, fin)

        for line in iterator:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = re.split(r"\s+", line.strip())
            if not _is_repeatmasker_data(parts):
                skipped.append(line.rstrip())
                continue
            try:
                query_start = int(parts[5]) - 1
                query_end = int(parts[6])
                consensus_start, consensus_end = _parse_consensus_interval(parts)
                perc_div = float(parts[1])
            except Exception:
                skipped.append(line.rstrip())
                continue
            if query_end <= query_start:
                skipped.append(line.rstrip())
                continue

            fragments.append(
                RepeatMaskerLocus(
                    query_name=parts[4],
                    query_start=query_start,
                    query_end=query_end,
                    family=parts[9],
                    repeat_class=parts[10],
                    strand="-" if parts[8] == "C" else "+",
                    identity=max(0.0, 100.0 - perc_div),
                    covered_bp=query_end - query_start,
                    consensus_start=consensus_start,
                    consensus_end=consensus_end,
                    rm_id=_parse_rm_id(parts),
                )
            )

    if log_path and skipped:
        with open(log_path, "w") as logf:
            logf.write("\n".join(skipped) + "\n")

    return fragments


def _merge_intervals(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
            continue
        total += cur_end - cur_start
        cur_start, cur_end = start, end
    total += cur_end - cur_start
    return total


@dataclass
class _ClusterBuilder:
    fragments: list[RepeatMaskerLocus]
    direction: str | None = None

    @property
    def query_start(self) -> int:
        return min(f.query_start for f in self.fragments)

    @property
    def query_end(self) -> int:
        return max(f.query_end for f in self.fragments)

    @property
    def consensus_start(self) -> int:
        return min(f.consensus_start for f in self.fragments)

    @property
    def consensus_end(self) -> int:
        return max(f.consensus_end for f in self.fragments)

    @property
    def last(self) -> RepeatMaskerLocus:
        return self.fragments[-1]


def _extends_right(
    builder: _ClusterBuilder,
    fragment: RepeatMaskerLocus,
    *,
    max_cons_gap: int,
    min_extension: int,
    max_backtrack: int,
) -> bool:
    return (
        fragment.consensus_end > builder.consensus_end + min_extension
        and fragment.consensus_start <= builder.consensus_end + max_cons_gap
        and fragment.consensus_start >= builder.consensus_start - max_backtrack
    )


def _extends_left(
    builder: _ClusterBuilder,
    fragment: RepeatMaskerLocus,
    *,
    max_cons_gap: int,
    min_extension: int,
    max_backtrack: int,
) -> bool:
    return (
        fragment.consensus_start < builder.consensus_start - min_extension
        and fragment.consensus_end >= builder.consensus_start - max_cons_gap
        and fragment.consensus_end <= builder.consensus_end + max_backtrack
    )


def _can_merge(
    builder: _ClusterBuilder,
    fragment: RepeatMaskerLocus,
    *,
    max_query_gap: int,
    max_cons_gap: int,
    min_extension: int,
    max_backtrack: int,
) -> tuple[bool, str | None]:
    if fragment.query_start - builder.query_end > max_query_gap:
        return False, None

    if (
        builder.last.rm_id is not None
        and fragment.rm_id is not None
        and builder.last.rm_id == fragment.rm_id
    ):
        if builder.direction == "+" and _extends_right(
            builder,
            fragment,
            max_cons_gap=max_cons_gap,
            min_extension=min_extension,
            max_backtrack=max_backtrack,
        ):
            return True, "+"
        if builder.direction == "-" and _extends_left(
            builder,
            fragment,
            max_cons_gap=max_cons_gap,
            min_extension=min_extension,
            max_backtrack=max_backtrack,
        ):
            return True, "-"
        if _extends_right(
            builder,
            fragment,
            max_cons_gap=max_cons_gap,
            min_extension=min_extension,
            max_backtrack=max_backtrack,
        ):
            return True, "+"
        if _extends_left(
            builder,
            fragment,
            max_cons_gap=max_cons_gap,
            min_extension=min_extension,
            max_backtrack=max_backtrack,
        ):
            return True, "-"
        return True, builder.direction

    if builder.direction in {None, "+"} and _extends_right(
        builder,
        fragment,
        max_cons_gap=max_cons_gap,
        min_extension=min_extension,
        max_backtrack=max_backtrack,
    ):
        return True, "+"
    if builder.direction in {None, "-"} and _extends_left(
        builder,
        fragment,
        max_cons_gap=max_cons_gap,
        min_extension=min_extension,
        max_backtrack=max_backtrack,
    ):
        return True, "-"
    return False, None


def _finalize_cluster(builder: _ClusterBuilder) -> RepeatMaskerLocus:
    fragments = list(builder.fragments)
    family_cov: dict[str, int] = defaultdict(int)
    family_ident: dict[str, float] = defaultdict(float)
    family_class: dict[str, str] = {}
    strand_cov = defaultdict(int)
    rm_ids = {f.rm_id for f in fragments if f.rm_id is not None}

    for fragment in fragments:
        span = fragment.length
        family_cov[fragment.family] += span
        family_ident[fragment.family] += fragment.identity * span
        family_class.setdefault(fragment.family, fragment.repeat_class)
        strand_cov[fragment.strand] += span

    best_family = max(
        family_cov,
        key=lambda family: (family_cov[family], family_ident[family]),
    )
    identity = family_ident[best_family] / family_cov[best_family]
    covered_bp = _merge_intervals([(f.query_start, f.query_end) for f in fragments])
    if builder.direction is None:
        if len(fragments) > 1:
            first = fragments[0]
            last = fragments[-1]
            if last.consensus_end > first.consensus_end:
                direction = "+"
            elif last.consensus_start < first.consensus_start:
                direction = "-"
            else:
                direction = "+" if strand_cov["+"] >= strand_cov["-"] else "-"
        else:
            direction = fragments[0].strand
    else:
        direction = builder.direction

    return RepeatMaskerLocus(
        query_name=fragments[0].query_name,
        query_start=min(f.query_start for f in fragments),
        query_end=max(f.query_end for f in fragments),
        family=best_family,
        repeat_class=family_class[best_family],
        strand=direction,
        identity=identity,
        covered_bp=covered_bp,
        consensus_start=min(f.consensus_start for f in fragments),
        consensus_end=max(f.consensus_end for f in fragments),
        rm_id=next(iter(rm_ids)) if len(rm_ids) == 1 else None,
        fragments=fragments,
    )


def collapse_repeatmasker_loci(
    fragments: list[RepeatMaskerLocus],
    *,
    max_query_gap: int = 250,
    max_cons_gap: int = 250,
    min_extension: int = 50,
    max_backtrack: int = 250,
) -> list[RepeatMaskerLocus]:
    """Collapse fragmented L1 RepeatMasker hits into locus-level calls."""

    by_query: dict[str, list[RepeatMaskerLocus]] = defaultdict(list)
    for fragment in fragments:
        by_query[fragment.query_name].append(fragment)

    collapsed: list[RepeatMaskerLocus] = []
    for query_name in sorted(by_query):
        hits = sorted(by_query[query_name], key=lambda hit: (hit.query_start, hit.query_end))
        current: _ClusterBuilder | None = None
        for hit in hits:
            if not hit.is_l1:
                if current is not None:
                    collapsed.append(_finalize_cluster(current))
                    current = None
                collapsed.append(hit)
                continue

            if current is None:
                current = _ClusterBuilder([hit])
                continue

            can_merge, direction = _can_merge(
                current,
                hit,
                max_query_gap=max_query_gap,
                max_cons_gap=max_cons_gap,
                min_extension=min_extension,
                max_backtrack=max_backtrack,
            )
            if can_merge:
                current.fragments.append(hit)
                if current.direction is None:
                    current.direction = direction
                continue

            collapsed.append(_finalize_cluster(current))
            current = _ClusterBuilder([hit])

        if current is not None:
            collapsed.append(_finalize_cluster(current))

    return sorted(collapsed, key=lambda hit: (hit.query_name, hit.query_start, hit.query_end))


def read_collapsed_repeatmasker_out(
    input_path: str | Path,
    *,
    log_path: str | Path | None = None,
    max_query_gap: int = 250,
    max_cons_gap: int = 250,
    min_extension: int = 50,
    max_backtrack: int = 250,
) -> list[RepeatMaskerLocus]:
    """Read a RepeatMasker ``.out`` file and collapse fragmented L1 loci."""

    fragments = read_repeatmasker_out(input_path, log_path=log_path)
    return collapse_repeatmasker_loci(
        fragments,
        max_query_gap=max_query_gap,
        max_cons_gap=max_cons_gap,
        min_extension=min_extension,
        max_backtrack=max_backtrack,
    )
