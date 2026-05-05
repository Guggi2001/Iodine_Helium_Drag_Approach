# Agent protocol

## Investigation mode

When the user asks to investigate, audit, inspect, compare, or explain:

- Do not edit files.
- Read relevant docs and tests.
- Report files inspected.
- Report conclusions and uncertainties.
- Suggest the smallest safe next edit.

## Edit mode

When the user explicitly asks to implement or fix:

- Make the smallest coherent change.
- Do not touch unrelated files.
- Add or update tests.
- Run relevant tests.
- Show changed files and remaining risks.

## Scientific-code caution

Clean code is not automatically correct physics.

Before changing a formula, unit conversion, force sign, random sampling method,
normalization convention, or indexing convention:

- locate the corresponding MATLAB source or previous Python test,
- explain the convention,
- add a regression test,
- only then edit.

## Forbidden without explicit user approval

- deleting reference data,
- changing physical constants,
- changing checkpoint schema,
- changing random-number draw order,
- changing default simulation scope,
- broad refactors,
- optimizing performance by changing numerical behavior,
- implementing out-of-scope MATLAB paths.