import pathlib
from datetime import datetime
from typing import Any, Dict, Union

from pydantic import AmqpDsn, BaseSettings, Field, MongoDsn, PostgresDsn, validator


class Settings(BaseSettings):
    """Settings"""

    APP_NAME: str = ""
    APP_VERSION: str = "1.0"
    APP_DESCRIPTION: str = ""

    LOG_LEVEL: str = Field(env="LOG_LEVEL", default="info")

    DB_PROVIDER: str = Field(env="DB_PROVIDER", default="postgres")
    DB_HOST: str = Field(env="DB_HOST", default="")
    DB_PORT: str = Field(env="DB_PORT", default="")
    DB_USER: str = Field(env="DB_USER", default="")
    DB_PASSWORD: str = Field(env="DB_PASSWORD", default="")
    DB_NAME: str = Field(env="DB_NAME", default="")
    DB_URI: Union[str, PostgresDsn, MongoDsn, None] = Field(
        env="DB_URI",
        default="postgresql://localhost/substances?user=admin&password=admin",
    )

    @classmethod
    @validator("DB_URI", pre=True)
    def assemble_db_connection(
        cls, value: Union[str, None], values: Dict[str, Any]
    ) -> str:
        """assemble_db_connection"""

        if isinstance(value, str):
            return value

        if values.get("DB_PROVIDER") == "postgres":
            return PostgresDsn.build(
                scheme="postgresql+asyncpg",
                user=values.get("DB_USER"),
                password=values.get("DB_PASSWORD"),
                host=values.get("DB_HOST"),
                port=values.get("DB_PORT"),
                path=f"/{0}".format(values.get("DB_NAME")),
            )

        raise ValueError("You must provide a valid DB_PROVIDER setting")

    RABBITMQ_HOST: str = Field(env="RABBITMQ_HOST", default="localhost")
    RABBITMQ_PORT: int = Field(env="RABBITMQ_PORT", default=5672)
    RABBITMQ_NAME: str = Field(env="RABBITMQ_NAME", default="substances")
    RABBITMQ_USER: str = Field(env="RABBITMQ_USER", default="")
    RABBITMQ_PASSWORD: str = Field(env="RABBITMQ_PASSWORD", default="")
    RABBITMQ_CONCURENT_CONSUMER: int = Field(
        env="RABBITMQ_CONCURENT_CONSUMER", default=1
    )
    RABBITMQ_URI: Union[AmqpDsn, None] = None

    @classmethod
    @validator("RABBITMQ_URI", pre=True)
    def assemble_rabbitmq_connection(
        cls,
        value: Union[str, None],
        values: Dict[str, Any],
    ) -> str:
        """assemble_rabbitmq_connection"""

        if isinstance(value, str):
            return value

        return AmqpDsn.build(
            scheme="amqp",
            user=values.get("RABBITMQ_USER"),
            password=values.get("RABBITMQ_PASSWORD"),
            host=values.get("RABBITMQ_HOST"),
            port=values.get("RABBITMQ_PORT"),
            path=f"/{0}".format(values.get("RABBITMQ_NAME")),
        )

    PUBCHEM_FTP_HOST: str = Field(
        env="PUBCHEM_FTP_HOST", default="ftp.ncbi.nlm.nih.gov"
    )
    PUBCHEM_FETCHER_ROOT_DIR: str = Field(env="PUBCHEM_FETCHER_ROOT_DIR", default=None)

    @classmethod
    @validator("PUBCHEM_FETCHER_ROOT_DIR", pre=True, always=True)
    def assemble_pubchem_root_dir(
        cls,
        value: Union[str, None],
        _: Dict[str, Any],
    ) -> str:
        if isinstance(value, str):
            return value

        now = datetime.now()
        return f'pubchem/Compound/Daily/{now.strftime("%Y-%m-%d")}'

    PUBCHEM_FETCHER_WORKING_DIR: str = Field(
        env="PUBCHEM_FETCHER_WORKING_DIR", default=None
    )

    @classmethod
    @validator("PUBCHEM_FETCHER_WORKING_DIR", pre=True, always=True)
    def assemble_pubchem_working_dir(
        cls,
        value: Union[str, None],
        _: Dict[str, Any],
    ) -> str:
        if isinstance(value, str):
            return value

        return str(pathlib.PurePosixPath(value))


settings = Settings()
