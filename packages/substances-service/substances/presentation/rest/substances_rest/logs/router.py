import asyncio
import json
import logging
from uuid import uuid4

from fastapi import Request, APIRouter
from sse_starlette import EventSourceResponse

router = APIRouter()
logger = logging.getLogger()
logs_queue = []


class LoggerInterceptor(logging.Handler):
    """LoggerInterceptor"""

    def emit(self, record: logging.LogRecord) -> None:
        logs_queue.append(record.getMessage())

log_interceptor = LoggerInterceptor()

@router.get("/stream")
async def log_stream(request: Request):
    """Retrieve all substances."""

    if log_interceptor not in logger.handlers:
        logger.addHandler(log_interceptor)

    async def event_generator():
        def dequeue_log():
            if len(logs_queue) > 0:
                return logs_queue.pop(0)
            return None

        while True:
            # If client closes connection, stop sending events
            if await request.is_disconnected():
                break

            try:
                log = dequeue_log()

                if log is not None:
                    yield {
                        "event": "log",
                        "id": str(uuid4()),
                        "data": json.dumps(
                            log,
                            default={"error": "JSON Serializion Failed"},
                            skipkeys=True,
                        ),
                    }

                await asyncio.sleep(0.3)
            except Exception as err:
                print(err)
                break

    return EventSourceResponse(event_generator())
