try:
    import os
    import sys

    # EFS - Add EFS mount path that contains libs into PYTHONPATH
    sys.path.append(os.environ["EFS_MOUNT_DIR"])  # nopep8 # noqa
except ImportError:
    pass

from substances_core.application.syncs.commands.execute import ExecuteSyncCommand

from .dependencies import execute_sync_command_handler  # type: ignore


def handle(event, context):  # type: ignore
    provider = event["provider"]
    data_location = event["data_location"] if "data_location" in event else None
    indexes = event["indexes"] if "indexes" in event else set()

    handler = execute_sync_command_handler()
    result = handler.execute(
        ExecuteSyncCommand(
            provider=provider, data_location=data_location, indexes=indexes
        )
    )

    return dict(result)
