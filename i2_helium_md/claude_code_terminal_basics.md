# Claude Code terminal basics

This guide is written for your current setup:

- Windows PC
- PyCharm terminal
- Project root opened at the folder containing `CLAUDE.md`
- Git repository already initialized
- Python migration project with `legacy_matlab_repository/`

The goal is not to memorize every Claude Code feature. The goal is to know the safe basic commands and how to use Claude Code without losing control of your codebase.

---

## 1. Where to start Claude Code

Always start Claude Code from the project root, not from a parent folder.

Good working directory:

```text
i2_helium_md/
```

This folder should contain:

```text
CLAUDE.md
README.md
current_state.md
next_tasks.md
migration_log.md
pyproject.toml
i2_helium_md/
tests/
docs/
legacy_matlab_repository/
```

Start Claude Code with:

```powershell
claude
```

If you are unsure where you are, run:

```powershell
pwd
dir
```

You should see `CLAUDE.md` in the directory listing.

---

## 2. Most important terminal commands outside Claude Code

These are normal terminal commands, not Claude Code prompts.

### Check Claude Code version

```powershell
claude --version
```

### Start Claude Code

```powershell
claude
```

### Start Claude Code with an initial prompt

```powershell
claude "explain this project structure"
```

### Continue the most recent conversation in the current directory

```powershell
claude -c
```

### One-shot question, then exit

```powershell
claude -p "explain what this repository does"
```

This is useful when you want a quick answer without entering a long interactive session.

### Check Git state before using Claude Code

```powershell
git status
```

### See changes after Claude Code edits files

```powershell
git diff
```

### Run tests

```powershell
pytest
```

Or run a specific test file:

```powershell
pytest tests/test_ion_propagation_step.py
```

---

## 3. Most important commands inside Claude Code

Inside an interactive Claude Code session, commands usually start with `/`.

Type:

```text
/
```

to see available commands in your installed version.

Availability can vary depending on version, platform, plan, plugins, and environment.

### Show help

```text
/help
```

Use this when you are unsure which command exists.

### Clear the current conversation context

```text
/clear
```

Use this when a session becomes confused or too long.

Important: before clearing, make sure important decisions are saved in files such as `current_state.md`, `next_tasks.md`, or `migration_log.md`.

### Compact the current conversation

```text
/compact
```

This summarizes older context to reduce context usage while keeping the session going.

Use this before the context gets too full.

### Check context/token usage

```text
/context
```

or check the context indicator shown in the Claude Code UI.

Use this when you suspect the session is getting large.

### Change or inspect model

```text
/model
```

Use this if you want to see or change the model used by Claude Code.

### Check permissions/settings

```text
/config
```

Use this to inspect or adjust Claude Code configuration and permission-related behavior.

### Exit Claude Code

Press:

```text
Ctrl+D
```

or type an exit command if your version shows one in `/help`.

---

## 4. Keyboard shortcuts worth knowing

### Cancel current generation or input

```text
Ctrl+C
```

Use this if Claude Code starts doing something you did not want.

### Exit session

```text
Ctrl+D
```

### See keyboard help

```text
?
```

Some shortcuts vary by terminal and operating system. The `?` key is the safest way to check what is available in your session.

---

## 5. Planning mode, edit approval, and permission modes

Claude Code has different permission modes. The most important beginner distinction is:

```text
Plan mode = Claude investigates and proposes a plan before editing.
Edit mode / default mode = Claude may ask for approval before edits.
Accept edits mode = Claude can apply file edits more automatically.
```

For your project, prefer this order:

```text
Plan first -> review plan -> approve one small implementation -> review diff -> run tests
```

### Switch modes during a session

In the Claude Code terminal, press:

```text
Shift+Tab
```

to cycle through modes.

The usual cycle is:

```text
default -> acceptEdits -> plan
```

The current mode is shown in the Claude Code status bar.

In some Windows terminal setups, the shortcut behavior can differ. If `Shift+Tab` does not show plan mode, press:

```text
?
```

inside Claude Code and check the keyboard shortcuts shown by your installed version.

You can also check:

```text
/help
```

or:

```text
/config
```

to inspect available commands and permission settings.

### What plan mode is for

Use plan mode when the task is non-trivial, for example:

- implementing the ion driver,
- designing MATLAB cross-reference tests,
- changing checkpoint handling,
- comparing Python and MATLAB behavior,
- refactoring a module,
- changing numerical logic,
- changing tests.

In plan mode, ask Claude Code to inspect files and propose a plan, not edit immediately.

Example:

```text
Plan only. Do not edit files.

Inspect the Python ion modules, the existing neutral driver, the relevant tests, and the MATLAB ion propagation source in legacy_matlab_repository/.

Propose a minimal implementation plan for Step 11d, the full ion propagation driver.

The plan must include:
1. files to inspect,
2. files likely to change,
3. tests to add,
4. risks,
5. MATLAB/Python deviations to watch for.
```

### How to approve a plan

After Claude Code proposes a plan, read it carefully.

If the plan is good, approve it explicitly with a narrow instruction:

```text
Approved. Implement only this plan.

Do not add extra features.
Do not change checkpoint schema.
Do not change physical constants.
Do not refactor unrelated files.
Run the relevant tests afterward.
```

If the plan is too broad, do not approve it. Instead say:

```text
Do not implement yet.

Reduce the plan to the smallest safe first step.
Only include the files needed for that step.
```

If you want to change the plan, say:

```text
Do not implement yet.

Modify the plan as follows:
- remove plotting,
- remove MATLAB cross-reference tests for now,
- only implement the ion driver and focused driver tests.
```

### How to reject or stop a plan

If Claude Code starts moving in the wrong direction, press:

```text
Ctrl+C
```

Then say:

```text
Stop. Do not edit more files.

Summarize what you changed so far and why.
```

Then inspect changes outside Claude Code:

```powershell
git diff
```

### Recommended mode for you

For scientific code, use plan mode often.

Recommended workflow:

```text
1. default mode for short explanations
2. plan mode for any implementation task
3. approve only one small implementation
4. avoid acceptEdits mode until you are confident
```

Do not use broad auto-approval behavior for physics code unless the change is very small and you already understand what it will do.

---

## 6. Safe first prompt for every new Claude Code session

Use this when starting a new session in this project:

```text
Read CLAUDE.md, README.md, current_state.md, next_tasks.md, testing.md, and migration_log.md.

Do not edit files.

Summarize:
1. the current migration state,
2. the next intended task,
3. which MATLAB legacy files are relevant,
4. which Python tests should be inspected first.

Do not modify files.
```

This verifies that Claude Code understands the project before touching anything.

---

## 7. Safe audit prompt

Use this before implementing anything:

```text
Audit the repository against the documentation.

Read CLAUDE.md, README.md, current_state.md, next_tasks.md, testing.md, and migration_log.md.

Inspect the relevant source files, but do not edit anything.

Report:
1. stale documentation,
2. mismatches between docs and code,
3. missing tests,
4. likely next safe task.

Do not modify files.
```

---

## 8. Safe implementation prompt pattern

Do not say:

```text
Finish the project.
```

Instead say:

```text
Implement only this task: <specific task>.

Before editing:
- inspect the relevant Python module,
- inspect the relevant Python tests,
- inspect the corresponding MATLAB source in legacy_matlab_repository/.

Constraints:
- do not change physical constants,
- do not change checkpoint schema,
- do not refactor unrelated files,
- do not implement out-of-scope features,
- add or update focused tests,
- run the relevant tests.

After editing, report:
1. files inspected,
2. MATLAB files inspected,
3. files changed,
4. tests run,
5. remaining risks.
```

---

## 9. Your recommended next coding prompt

For your current project state, the next coding task should likely be the full ion driver.

Use this only after an audit looks reasonable:

```text
Implement Step 11d only: the full ion propagation driver.

Follow CLAUDE.md and next_tasks.md.

Before editing, inspect:
- i2_helium_md/simulation/neutral.py
- i2_helium_md/simulation/ion_initial_state.py
- i2_helium_md/simulation/ion_propagation_step.py
- i2_helium_md/simulation/checkpoint.py
- i2_helium_md/simulation/run_directory.py
- the relevant tests
- the corresponding MATLAB ion propagation source in legacy_matlab_repository/

Do not implement MATLAB cross-reference tests yet.
Do not implement plotting.
Do not implement the single-pulse script yet.
Do not change checkpoint schema.
Do not change physical constants.
Do not refactor neutral propagation unless required by a failing test.

Add focused tests for the ion driver and run the relevant tests.
```

---

## 10. How to reference files in prompts

You can explicitly mention files by path:

```text
Inspect i2_helium_md/simulation/neutral.py and compare it with i2_helium_md/simulation/ion_propagation_step.py.
```

For your MATLAB reference:

```text
Search legacy_matlab_repository/ for vmi_sim_3d_ion_propa.m and inspect the ion propagation loop.
```

Keep paths simple with forward slashes:

```text
i2_helium_md/simulation/ion.py
legacy_matlab_repository/
testing.md
```

Avoid Windows backslashes in prompts unless necessary.

---

## 11. When to use /clear vs /compact

Use `/compact` when:

- the session is useful,
- Claude Code understands the task,
- context is getting large,
- you want to continue without starting over.

Use `/clear` when:

- Claude Code seems confused,
- the task changed completely,
- the session is polluted by old assumptions,
- you want a fresh start.

Before `/clear`, update project files if important decisions were made.

---

## 12. Beginner safety workflow

Use this pattern:

```text
1. git status
2. claude
3. switch to plan mode for implementation tasks
4. audit-only prompt
5. ask for implementation plan
6. approve one small implementation
7. let Claude Code edit
8. git diff
9. run tests
10. manually review changes
11. commit yourself only when you understand the diff
```

Do not let Claude Code make a big multi-feature change in one step.

---

## 13. Commands you should avoid as a beginner

Avoid starting Claude Code with dangerous broad permission modes unless you know exactly what they do.

Do not use commands or flags that skip permission prompts unless explicitly needed.

Avoid prompts like:

```text
Delete unused files.
Clean up the whole project.
Refactor everything.
Make the project better.
Fix all issues.
Finish the migration.
```

These are too broad for scientific code.

---

## 14. What to do if Claude Code starts editing too much

Press:

```text
Ctrl+C
```

Then say:

```text
Stop. Do not edit more files.

Summarize exactly what you changed so far and why.
```

Then inspect:

```powershell
git diff
```

If the changes are bad and you want to discard all uncommitted changes:

```powershell
git restore .
```

Be careful: `git restore .` deletes all uncommitted edits in tracked files.

If there are new untracked files and you want to see them:

```powershell
git status
```

Do not delete untracked files automatically unless you are sure.

---

## 15. Minimal command cheat sheet

### Outside Claude Code

```powershell
claude
claude --version
claude -c
claude -p "question"
git status
git diff
pytest
pytest tests/test_file.py
```

### Inside Claude Code

```text
/
/help
/clear
/compact
/context
/model
/config
Shift+Tab
Ctrl+C
Ctrl+D
?
```

---

## 16. Best mindset

Claude Code is not a magic senior scientist.

Use it like a careful junior developer:

- give narrow tasks,
- use plan mode before implementation,
- require tests,
- require MATLAB cross-reference,
- review every diff,
- do not accept physics changes you do not understand.

For this project, the safest sequence is:

```text
audit -> ion driver -> ion MATLAB cross-reference tests -> run script -> HeDFT loader -> trajectory comparison
```
