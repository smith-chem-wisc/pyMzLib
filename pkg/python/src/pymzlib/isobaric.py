"""Isobaric labelling kits - TMT, TMTpro, iTRAQ, DiLeu - with every channel's reporter-ion m/z.

An isobaric experiment is read out of the low-m/z reporter ions: one ion per channel, a few
millidaltons apart. Getting a channel's label or its m/z wrong mislabels every sample after it, so
pyMzLib does not keep its own table. :func:`kits` returns mzLib's::

    >>> import pymzlib
    >>> tmtpro = pymzlib.isobaric.kits("TMT18")
    >>> kit = tmtpro.kits[0]
    >>> kit.name, kit.channel_count
    ('TMT18', 18)
    >>> for channel in kit.channels[:3]:
    ...     print(channel.label, round(channel.reporter_ion_mz, 5))
    126 126.12773
    127N 127.12476
    127C 127.13108

**Nothing in that table is typed in.** mzLib holds only the channel labels. Each m/z is a ``DI HCD``
diagnostic-ion line of the kit's ``Multiplex Label`` modification in mzLib's embedded ``TMT.txt``,
plus one proton, sorted ascending and paired with the labels by position. The m/z are therefore
the same ones MetaMorpheus uses to quantify a TMT search.

The kit is looked up the way MetaMorpheus looks it up: by its whole name, case-insensitively,
optionally followed by ``" on <motif>"`` - ``"TMT10"``, ``"tmt10"``, ``"TMT6-plex"``,
``"iTRAQ-4plex on K"``. Never by substring, so ``"TMT10plex"`` is refused rather than mistaken for a
kit whose name it contains.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

from . import _bridge
from .readers import _Table

__all__ = [
    "ReporterChannel",
    "IsobaricKit",
    "IsobaricKits",
    "kits",
]

@dataclass(frozen=True)
class ReporterChannel:
    """One channel of a kit and the theoretical m/z of its reporter ion.

    Attributes:
        kit: mzLib's kit name (an ``IsobaricMassTagType`` member, e.g. ``"TMT18"``).
        index: 0-based position within the kit, ascending with ``reporter_ion_mz``.
        label: The channel's name as MetaMorpheus spells it - ``"126"``, ``"127N"``, ``"115a"``.
        reporter_ion_mz: Theoretical m/z at charge 1, in **m/z** (Th): the diagnostic-ion mass
            from mzLib's ``TMT.txt`` plus one proton. Not a calibrated or observed value.
        mz_min: ``reporter_ion_mz`` minus the matching tolerance, in m/z.
        mz_max: ``reporter_ion_mz`` plus the matching tolerance, in m/z.
    """

    kit: str
    index: int
    label: str
    reporter_ion_mz: float
    mz_min: float
    mz_max: float


@dataclass(frozen=True)
class IsobaricKit:
    """One kit and its channels, in ascending reporter-ion m/z.

    Attributes:
        name: mzLib's kit name, an ``IsobaricMassTagType`` member. ``TMT16`` and ``TMT18`` are
            TMTpro; ``kits().kits`` lists every name.
        channel_count: How many channels the kit has.
        channels: Every channel as a :class:`ReporterChannel`.
    """

    name: str
    channel_count: int
    channels: list[ReporterChannel]

    @property
    def labels(self) -> list[str]:
        """The channel labels, in ascending reporter-ion m/z."""
        return [c.label for c in self.channels]

    @property
    def reporter_ion_mzs(self) -> list[float]:
        """The reporter-ion m/z values, ascending, in m/z."""
        return [c.reporter_ion_mz for c in self.channels]


@dataclass(frozen=True)
class IsobaricKits(_Table):
    """The kits :func:`kits` returned, as a long table and as objects.

    ``columns`` maps ``kit``, ``channel_index``, ``channel_label``, ``reporter_ion_mz``, ``mz_min``
    and ``mz_max`` to one value per (kit, channel) - ready for ``pandas.DataFrame(result.columns)``.
    :attr:`kits` groups the same rows by kit.

    Attributes:
        kit: The name you asked for, exactly as given, or ``None`` when every kit was listed.
        kit_count: Kits listed.
        record_count: Channels listed, over every kit; the length of every column.
        absolute_tolerance: The half-width, in **Da**, of each channel's matching window. mzLib
            reads a reporter intensity as the most intense peak within it.
        column_names: The table's columns, in order.
        columns: Column name to one value per channel.
        kit_summaries: The wire's per-kit summaries, ``{"kit", "channel_count"}``, in mzLib's
            order. :attr:`kits` regroups the table under them as :class:`IsobaricKit` objects.
        caveats: What the m/z are and are not - read these once.
    """

    kit: str | None
    kit_count: int
    record_count: int
    absolute_tolerance: float
    column_names: list[str]
    columns: dict[str, list[Any]]
    kit_summaries: list[dict[str, Any]] = field(default_factory=list, repr=False)
    caveats: list[str] = field(default_factory=list)

    @property
    def channels(self) -> list[ReporterChannel]:
        """Every channel as a :class:`ReporterChannel`, kits in mzLib's order."""
        return [
            ReporterChannel(
                kit=row["kit"],
                index=row["channel_index"],
                label=row["channel_label"],
                reporter_ion_mz=row["reporter_ion_mz"],
                mz_min=row["mz_min"],
                mz_max=row["mz_max"],
            )
            for row in self.records
        ]

    @property
    def kits(self) -> list[IsobaricKit]:
        """Every kit as an :class:`IsobaricKit`, in mzLib's order."""
        channels = self.channels
        return [
            IsobaricKit(
                name=summary["kit"],
                channel_count=summary["channel_count"],
                channels=[c for c in channels if c.kit == summary["kit"]],
            )
            for summary in self.kit_summaries
        ]

    def kit_named(self, name: str) -> IsobaricKit:
        """The listed kit whose :attr:`IsobaricKit.name` is exactly ``name``.

        Raises:
            UsageError: no listed kit has that name.
        """
        for kit in self.kits:
            if kit.name == name:
                return kit
        raise _bridge.UsageError(
            f"No kit named {name!r} in this result; it lists {[k.name for k in self.kits]}."
        )

    @classmethod
    def _from_wire(cls, data: Mapping[str, Any]) -> "IsobaricKits":
        return cls(
            kit=data.get("kit"),
            kit_count=int(data.get("kit_count", 0)),
            record_count=int(data.get("record_count", 0)),
            absolute_tolerance=float(data.get("absolute_tolerance", 0.0)),
            column_names=list(data.get("column_names") or []),
            columns={k: list(v) for k, v in (data.get("columns") or {}).items()},
            kit_summaries=[dict(k) for k in (data.get("kits") or [])],
            caveats=list(data.get("caveats") or []),
        )


def kits(kit: str | None = None, *, timeout: float | None = 60) -> IsobaricKits:
    """The isobaric kits mzLib can name, with every channel's label and reporter-ion m/z.

    Calls mzLib's ``IsobaricMassTag.TryGetIsobaricMassTag`` for each ``IsobaricMassTagType``, or for
    the one ``kit`` resolves to through ``IsobaricMassTag.TryGetTagType``. No file, no network.

    Args:
        kit: One kit, matched by mzLib's whole-name rule: case-insensitive, optionally followed by
            ``" on <motif>"`` (``"TMT10"``, ``"TMT6-plex"``, ``"iTRAQ-4plex on K"``). ``None`` lists
            every kit.
        timeout: Seconds to allow. ``None`` waits indefinitely.

    Returns:
        An :class:`IsobaricKits` table, also grouped as :attr:`IsobaricKits.kits`.

    Raises:
        UsageError: ``kit`` is blank or not a name mzLib knows. The message lists the kits.

    Examples:
        Every kit mzLib knows:

        >>> catalogue = kits()
        >>> [(k.name, k.channel_count) for k in catalogue.kits]  # doctest: +NORMALIZE_WHITESPACE
        [('TMT6', 6), ('TMT10', 10), ('TMT11', 11), ('TMT16', 16), ('TMT18', 18),
         ('iTRAQ4', 4), ('iTRAQ8', 8), ('diLeu4', 4), ('diLeu12', 12)]
        >>> catalogue.absolute_tolerance
        0.003
        >>> catalogue.kit_named("iTRAQ8").labels
        ['113', '114', '115', '116', '117', '118', '119', '121']

        TMTpro 18-plex, whose top channel sits at 135.15:

        >>> tmtpro = kits("TMT18")
        >>> last = tmtpro.kits[0].channels[-1]
        >>> last.label, round(last.reporter_ion_mz, 4), round(last.mz_max - last.mz_min, 6)
        ('135N', 135.1516, 0.006)
    """
    args = ["isobaric", "kits"]
    if kit is not None:
        if not isinstance(kit, str) or not kit.strip():
            raise _bridge.UsageError(
                f"kit must be a kit name such as 'TMT10' or 'iTRAQ-4plex', or None; got {kit!r}."
            )
        args += ["--kit", kit]
    data = _bridge.invoke(*args, timeout=timeout)
    return IsobaricKits._from_wire(data)
