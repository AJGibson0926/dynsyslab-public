"""Source-checkout entry point: ``python -m dynsyslab`` runs the front door.

The installable Python package is ``dynsyslab_client``, a thin client over the compiled suite.
"""

if __name__ == "__main__":
    from operations.front_door import main
    raise SystemExit(main())
