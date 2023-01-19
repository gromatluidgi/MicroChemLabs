from fastapi import APIRouter, Depends
from shared.application.commands import CommandHandler
from substances_core.application.substances.commands.import_string import (
    ImportSubstanceCommand,
    ImportSubstanceResult,
)

from .dependencies import import_substance_from_string_command_handler

router = APIRouter()


@router.get(
    "",
)
def get_substances():
    """Retrieve all substances."""
    return []


@router.get("/{inchikey}/render")
def render(inchikey: str):
    return "render"


@router.post(
    "/import",
)
def import_from_string(
    command: ImportSubstanceCommand,
    handler: CommandHandler = Depends(import_substance_from_string_command_handler),
) -> ImportSubstanceResult:
    return handler.execute(command)
