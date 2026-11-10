"""結果サイトを作ってローカルで公開する（./tasks.ps1 server と同じ）。

  python serve.py              # サイトを作って http://127.0.0.1:8000/ で公開、ブラウザで開く
  python serve.py --no-build   # 作り直さずに公開だけ
  python serve.py --port 8080
"""
from __future__ import annotations

import argparse
import functools
import http.server
import subprocess
import sys
import webbrowser

from common.runner import ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--no-open", action="store_true")
    ap.add_argument("--n3d", type=int, default=30)
    args = ap.parse_args()

    if not args.no_build:
        subprocess.run([sys.executable, str(ROOT / "build_site.py"), "--n3d", str(args.n3d)], check=True)

    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT / "results"))
    url = f"http://127.0.0.1:{args.port}/"
    with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler) as srv:
        print(f"serving {ROOT / 'results'} at {url}  (Ctrl+C で終了)")
        if not args.no_open:
            webbrowser.open(url)
        try:
            srv.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
