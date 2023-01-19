from substances_core.application.substances.commands.import_string import (
    ImportSubstanceCommand,
)

from .dependencies import import_substance_from_string_command_handler


def handle(event, context):
    print(event["data"])

    handler = import_substance_from_string_command_handler()
    result = handler.execute(
        ImportSubstanceCommand(data=event["data"], format=event["format"])
    )
    response = {"success": True, "body": result.json()}

    return response
