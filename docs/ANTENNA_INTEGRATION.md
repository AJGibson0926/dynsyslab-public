# Shared Antenna integration

dynsyslab uses the centralized **Antenna** repository for remote operations:
[`AJGibson0926/the-antenna`](https://github.com/AJGibson0926/the-antenna).
It does not vendor, fork, copy, modify, or submodule the Antenna source.

## Link resolution

`operations.antenna_link.AntennaLink` resolves the shared checkout from
`ANTENNA_HOME`, verifies its exact commit against
`config/antenna_dependency.json`, then executes the upstream launcher
`<ANTENNA_HOME>/tools/antenna.py`.

The shared checkout is currently installed at
`~/.codex/tools/the-antenna`, and that location has been set as
the user-level `ANTENNA_HOME`. A new shell picks up that setting automatically.

## dynsyslab policy

`config/antenna_policy.json` is owned by dynsyslab and versioned here. It
declares only the dynsyslab Git repository. The Antenna itself remains
default-deny: a resource and each allowed action must be named by policy.

- `inspect` is safe, read-only repository status.
- `sync` pushes only a clean, already-committed branch.
- `commit-push` requires exact paths, `--apply`, and the review rule recorded
  in the policy.

No S3 resource is active yet. `config/antenna_policy.s3.template.json` is a
non-operational starting point for when a bucket and AWS profile are supplied.
Never put credentials in either policy.

## Use

```powershell
# Calls the centrally managed Antenna; it does not use a dynsyslab copy.
python -m operations.antenna_link doctor
python -m operations.antenna_link policy-validate --policy config/antenna_policy.json
python -m operations.antenna_link git inspect --policy config/antenna_policy.json --resource dynsyslab-source
```

The wrapper is a link only. It deliberately has no remote-storage, compute, or
Git implementation of its own.
