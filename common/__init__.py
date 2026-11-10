import sys

# Windows のコンソール (cp932) で表示できない文字（µ など）があっても止まらないようにする
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except (AttributeError, ValueError):
        pass
