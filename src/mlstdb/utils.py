import click

# Color formatting using click
def error(msg: str) -> None:
    click.secho(f"Error: {msg}", fg="red", err=True)

def success(msg: str) -> None:
    click.secho(msg, fg="green")

def info(msg: str) -> None:
    click.secho(msg, fg="blue")


