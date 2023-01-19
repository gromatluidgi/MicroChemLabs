import logging
from typing import List, Optional

from fastapi import APIRouter
from pydantic import BaseModel
from starlette.responses import JSONResponse

from .logs.router import router as logs_router
from .substances.router import router as substances_router
from .syncs.router import router as syncs_router

LOGGER = logging.getLogger(__name__)


class ErrorMessage(BaseModel):
    """ErrorMessage"""

    msg: str


class ErrorResponse(BaseModel):
    """ErrorResponse"""

    detail: Optional[List[ErrorMessage]]


api_router = APIRouter(
    default_response_class=JSONResponse,
    responses={
        400: {"model": ErrorResponse},
        401: {"model": ErrorResponse},
        403: {"model": ErrorResponse},
        404: {"model": ErrorResponse},
        500: {"model": ErrorResponse},
    },
)

api_router.include_router(logs_router, prefix="/logs", tags=["logs"])
api_router.include_router(substances_router, prefix="/substances", tags=["substances"])
api_router.include_router(syncs_router, prefix="/syncs", tags=["syncs"])


@api_router.get("/healthcheck", include_in_schema=False)
def healthcheck():
    """healthcheck"""
    LOGGER.info("API Call /healthcheck")
    return {"status": "ok"}
