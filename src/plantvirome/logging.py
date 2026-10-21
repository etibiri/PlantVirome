"""JSON-lines logging to stderr and, optionally, a file."""
import json
import logging
import sys
import time


class JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        entry = {
            "ts": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(record.created)),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        entry.update(getattr(record, "fields", {}))
        return json.dumps(entry, default=str)


def setup_logging(logfile=None, level: int = logging.INFO) -> logging.Logger:
    log = logging.getLogger("plantvirome")
    log.setLevel(level)
    log.handlers.clear()
    for handler in (logging.StreamHandler(sys.stderr),
                    *([logging.FileHandler(logfile)] if logfile else [])):
        handler.setFormatter(JsonFormatter())
        log.addHandler(handler)
    return log


def event(msg: str, **fields) -> None:
    """Log a message with structured fields."""
    logging.getLogger("plantvirome").info(msg, extra={"fields": fields})
