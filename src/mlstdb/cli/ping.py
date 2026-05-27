import json

import click

from mlstdb.core.auth import ping_url


@click.command()
@click.help_option("-h", "--help")
@click.argument("url")
@click.option(
    "--db",
    "-d",
    type=click.Choice(["pubmlst", "pasteur"]),
    default=None,
    help="Database whose stored credentials are used for authentication.",
)
@click.option(
    "--no-auth",
    "no_auth",
    is_flag=True,
    default=False,
    help="Skip authentication and send an unauthenticated request.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help="Print extra request detail (auth mode, full URL, status code).",
)
def ping(url: str, db: str, no_auth: bool, verbose: bool):
    """Probe an API endpoint and display the response.

    \b
    Makes a single GET request to URL and prints the HTTP status and
    response body.  Authentication is attempted automatically in this order:

    \b
      1. API key  (stored by 'mlstdb connect --api-key')
      2. OAuth session token  (stored by 'mlstdb connect')
      3. Unauthenticated fallback  (with a warning)

    Use --no-auth to skip all authentication, equivalent to the
    --no-auth flag on 'mlstdb update'.

    \b
    Examples:
      mlstdb ping https://rest.pubmlst.org/db --no-auth
      mlstdb ping https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes --db pubmlst
      mlstdb ping https://bigsdb.pasteur.fr/api/db --db pasteur --verbose
    """
    # Prompt for --db when authentication is needed and it was not supplied
    if not no_auth and db is None:
        db = click.prompt(
            "Which database holds your credentials for this URL?",
            type=click.Choice(["pubmlst", "pasteur"]),
            default="pubmlst",
        )

    try:
        status_code, body, json_payload = ping_url(
            url=url,
            db=db,
            verbose=verbose,
            no_auth=no_auth,
        )
    except Exception as exc:
        error(f"\n{exc}")
        raise SystemExit(1)

    click.echo(f"\nHTTP {status_code}  {url}\n")

    if json_payload is not None:
        click.echo(json.dumps(json_payload, indent=2))
    else:
        click.echo(body)

    if status_code in (401, 403):
        click.echo(
            "\nHint: run 'mlstdb connect --db <db>' to configure credentials, "
            "or use --no-auth for unauthenticated access."
        )
