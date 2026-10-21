"""Run external tools with logged command, return code and duration."""
import shutil
import subprocess
import time
from pathlib import Path

from .logging import event


class ToolError(RuntimeError):
    pass


def require(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise ToolError(f"missing tools on PATH: {', '.join(missing)}")


def run(cmd: list[str], *, stdout: Path | None = None, cwd: Path | None = None) -> None:
    """Run `cmd`; raise ToolError on non-zero exit. `stdout` redirects to a file."""
    require(cmd[0])
    event("run", cmd=" ".join(map(str, cmd)))
    start = time.monotonic()
    out = open(stdout, "w") if stdout else None
    try:
        proc = subprocess.run([str(c) for c in cmd], stdout=out, stderr=subprocess.PIPE,
                              text=True, cwd=cwd)
    finally:
        if out:
            out.close()
    event("done", cmd=str(cmd[0]), rc=proc.returncode,
          seconds=round(time.monotonic() - start, 1))
    if proc.returncode:
        raise ToolError(f"{cmd[0]} failed (rc={proc.returncode}):\n{proc.stderr[-2000:]}")


def step_done(outdir: Path) -> bool:
    return (outdir / ".done").exists()


def mark_done(outdir: Path) -> None:
    (outdir / ".done").touch()
