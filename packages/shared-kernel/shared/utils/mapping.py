from dataclasses import _MISSING_TYPE, fields
from typing import Any, Optional

import pydantic


class MappingUtils:
    # Taken from: https://towardsdatascience.com/pydantic-or-dataclasses-why-not-both-convert-between-them-ba382f0f9a9c
    @staticmethod
    def convert_flat_dataclass_to_pydantic(
        dcls: type, name: Optional[str] = None
    ) -> type[pydantic.BaseModel]:
        if name is None:
            name_ = f"Pydantic{dcls.__name__}"
        else:
            name_ = name

        # type: ignore
        return pydantic.create_model(
            name_,
            **MappingUtils._get_pydantic_field_kwargs(dcls),
        )

    @staticmethod
    def _get_pydantic_field_kwargs(dcls: type) -> dict[str, tuple[type, Any]]:
        # get attribute names and types from dataclass into pydantic format
        pydantic_field_kwargs = dict()
        for _field in fields(dcls):
            # check is field has default value
            if isinstance(_field.default, _MISSING_TYPE):
                # no default
                default = ...
            else:
                default = _field.default

            pydantic_field_kwargs[_field.name] = (_field.type, default)
        return pydantic_field_kwargs
