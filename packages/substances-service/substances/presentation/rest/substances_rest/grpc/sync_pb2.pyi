from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Iterable as _Iterable, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class GetRequest(_message.Message):
    __slots__ = ["id"]
    ID_FIELD_NUMBER: _ClassVar[int]
    id: str
    def __init__(self, id: _Optional[str] = ...) -> None: ...

class ItemResponse(_message.Message):
    __slots__ = ["id", "provider"]
    ID_FIELD_NUMBER: _ClassVar[int]
    PROVIDER_FIELD_NUMBER: _ClassVar[int]
    id: str
    provider: str
    def __init__(self, id: _Optional[str] = ..., provider: _Optional[str] = ...) -> None: ...

class ListRequest(_message.Message):
    __slots__ = ["limit", "offset", "ordering"]
    LIMIT_FIELD_NUMBER: _ClassVar[int]
    OFFSET_FIELD_NUMBER: _ClassVar[int]
    ORDERING_FIELD_NUMBER: _ClassVar[int]
    limit: int
    offset: int
    ordering: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, offset: _Optional[int] = ..., limit: _Optional[int] = ..., ordering: _Optional[_Iterable[str]] = ...) -> None: ...

class ListResponse(_message.Message):
    __slots__ = ["count"]
    COUNT_FIELD_NUMBER: _ClassVar[int]
    count: int
    def __init__(self, count: _Optional[int] = ...) -> None: ...
