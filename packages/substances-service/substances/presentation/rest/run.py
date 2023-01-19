import uvicorn


def start():
    """Start application."""
    uvicorn.run(
        app="substances_rest.main:api",
        host="0.0.0.0",
        port=8080,
        reload=True,
        workers=1,
    )


if __name__ == "__main__":
    start()
