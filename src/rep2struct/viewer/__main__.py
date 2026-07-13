import sys

from . import serve


def main():
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print("usage: python -m rep2struct.viewer RUN_DIR [--port N]\n\n"
              "Serves a live timeline of RUN_DIR/transcript.jsonl in the browser.")
        return
    port = 8000
    if "--port" in args:
        i = args.index("--port")
        port = int(args[i + 1])
        del args[i:i + 2]
    serve(args[0], port=port)


if __name__ == "__main__":
    main()
