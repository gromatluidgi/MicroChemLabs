from substances_core.application.syncs.commands.execute import ExecuteSyncCommand


def test_handle():
    # Arrange
    provider = "PUBCHEM"
    command = ExecuteSyncCommand(provider=provider)

    # Act

    # Assert
