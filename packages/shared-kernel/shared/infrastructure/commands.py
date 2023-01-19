from typing import Callable, Dict

from shared.application.commands import (
    CommandBase,
    CommandDispatcher,
    CommandHandler,
    CommandResultBase,
)


class InternalCommandDispatcher(CommandDispatcher):
    """InternalCommandDispatcher"""

    def __init__(self) -> None:
        self._handlers: Dict[type[CommandBase], Callable] = {}

    def register_handler(
        self,
        command: type[CommandBase],
        handler: CommandHandler[CommandBase, CommandResultBase],
    ) -> None:
        if command not in self._handlers:
            self._handlers[command] = handler

    def dispatch(self, command: CommandBase):
        if type(command) in self._handlers:
            self._handlers[type(command)](command)
