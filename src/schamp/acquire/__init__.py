"""`schamp-acquire`: the run sheet an operator keeps beside the instrument.

`model` is the run sheet itself -- the planned series, the record it writes, the checks
it runs -- and imports no Qt, so everything it does is checkable without a display.
`app` and the widget modules beside it are the only places that import PySide6, which
is the optional extra `schamp[acquire]`; importing this package does not need it.
"""
