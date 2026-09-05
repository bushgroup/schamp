# Test fixtures

Small files the self-check needs. Nothing here is generated, and nothing here is large.

## `_extern.inf`

A verbatim copy of the `_extern.inf` sidecar of one 2013-04-23 poly-DL-alanine
acquisition on the UW Synapt G2: 192 lines, 5.4 kB, the instrument settings MassLynx
wrote when the file was acquired. It is the fixture for `schamp.extern`, and it is a
real one on purpose, because every interesting property of that parser is a property of
what the instrument actually writes rather than of the format anyone would design:

- **CRLF line endings and Latin-1 bytes.** The degree sign in `Source Temperature (°C)`
  is a raw `0xb0`, which is not valid UTF-8, so the parser's encoding fallback is
  exercised only if the bytes survive. The repo's `.gitattributes` pins `* text=auto
  eol=lf`; this one path is exempted with `-text` so the file stays byte for byte what
  MassLynx produced.
- **Keys padded out with tabs**, and section headers that are bare lines with no tab.
- **No `Pusher Interval` key.** `Pusher` is the pusher pulse voltage and `Pusher Offset`
  is a voltage too. The 2013 analysis grepped for a pusher interval, silently defaulted
  it to 1.0 and so produced every centroid in drift bins rather than milliseconds; the
  check asserts the key's absence so that nothing ever reintroduces that fallback.
- **Real drift-voltage optics**, so the instrument profile's formula can be evaluated
  against an acquisition rather than against a hand-written dictionary.

This is acquisition parameters and nothing else: no sample identities, no results, no
part of the Waters SDK. A `.raw` acquisition itself is hundreds of megabytes and is
never in this repository.

## `_HEADER.TXT`

**Synthetic, unlike its neighbour, and deliberately so.** Eleven lines in the layout
MassLynx writes -- `$$ Key: value`, CRLF -- carrying what the parser and
`schamp.calibration` have to handle and what no single real file holds all of at once:

- **`Cal Function 1`, a real cubic**, copied verbatim from a 2013-04-23 acquisition on
  the UW Synapt G2. It is the calibration polynomial the spectrum readers apply and the
  chromatogram reader does not, so it is the one line in this file that no SDK accessor
  exposes and the reason the parser exists at all.
- **`Cal Function 2`, the identity** `0,1`, which is what MassLynx writes for a function
  with no calibration of its own and which must stay distinguishable from an absent
  line.
- **`Cal Function 3`, unparseable on purpose**, so the check can assert that a
  half-readable calibration comes back as nothing rather than as a number.
- **`Cal StdDev Function 1`**, whose key begins the same way as the ones above, so a
  prefix match instead of an exact one would be caught.
- **A value containing colons** (`Cal Time: 13:36`) and **a line with no `$$`**, the two
  ways a naive split gets this format wrong.

Nothing else in it is real: no sample, no acquisition name, no operator. That is the
whole reason it is synthetic -- a real header carries sample descriptions and acquired
names, and this file needs none of them to test what it tests. The behaviour on real
files is pinned lab-side, against every acquisition of the 2013 series.
