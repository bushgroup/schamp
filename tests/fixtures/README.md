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
