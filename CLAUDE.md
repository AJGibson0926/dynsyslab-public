# DynSysLab for Claude Code

The engineering policy is shared with every agent and lives in `AGENTS.md`;
it is canonical. Do not fork it here.

@AGENTS.md

## First session

Follow "First session with an agent" in `AGENTS.md`: if `.dynsyslab/agent-welcomed` is
missing, greet the user with the overview in `docs/AGENT_WELCOME.md`, then create the file.

## Orientation

Read, in order, before editing: `docs/START_HERE.md`, `ARCHITECTURE.md`,
`docs/AI_CONTEXT.md`, `docs/AGENT_HARNESS.md`, then the relevant feature guide
under `docs/` or `docs/architecture/`. Use the `dynsyslab-spectral-analysis`
skill for DMD/Koopman/spectral KROM work.

The `dynsyslab-harness` MCP server (`codex_harness/server.py`) is agent-neutral
despite its name; it is read-only orientation only.

## Working alongside Codex (author's setup)

This section describes the author's own machine; in a public clone there is no Hub and no
other agent, and it can be ignored.

Codex also works in this checkout. Its front end (`AGENTS.md`, the routing Hub
at `Documents\Codex\DynSysLab-Project-Hub`, and `~/.codex`) must keep working.

- Before editing, run `git status`. Uncommitted changes you did not make are
  another agent's work in progress: do not modify, revert, stash, or commit
  them. Stop and ask the user.
- Do not edit the Hub or `~/.codex`. Change `AGENTS.md` only when the user asks
  for a policy change, since it governs both agents.
- Keep commits focused on your own task so either agent can review history.
