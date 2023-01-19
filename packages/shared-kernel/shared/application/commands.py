from abc import ABC, abstractmethod
from typing import Generic, TypeVar

from pydantic import BaseModel, validator


class CommandBase(BaseModel):
    """CommandBase"""

    command_name: str = ""

    @validator("command_name", pre=True, always=True)
    def class_name(cls, value):  # pylint: disable=no-self-argument
        if value is not None:
            return value
        return cls.__name__  # pylint: disable=no-member


class CommandResultBase(BaseModel):
    """CommandResultBase"""


# pylint: disable=invalid-name
TCommandBase = TypeVar("TCommandBase", bound=CommandBase)
TCommandResultBase = TypeVar("TCommandResultBase", bound=CommandResultBase)


class CommandHandler(ABC, Generic[TCommandBase, TCommandResultBase]):
    """CommandHandler"""

    @abstractmethod
    def execute(self, command: TCommandBase) -> TCommandResultBase:
        """Handle"""
        raise NotImplementedError

    def __call__(self, command: TCommandBase) -> TCommandResultBase:
        """Make this class a callable."""
        return self.execute(command)


class CommandDispatcher(ABC, Generic[TCommandBase]):
    """CommandDispatcher"""

    @abstractmethod
    def register_handler(
        self,
        command: type[TCommandBase],
        handler: CommandHandler[TCommandBase, TCommandResultBase],
    ) -> None:
        raise NotImplementedError

    @abstractmethod
    def dispatch(self, command: TCommandBase) -> None:
        raise NotImplementedError
