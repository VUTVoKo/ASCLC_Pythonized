# Project instructions

## Communication rules

- Only write or change code when explicitly told to. By default, requests to
  look at, explain, investigate, or discuss something call for an answer, not
  an edit. Wait for an instruction to implement.
- When the user tells you something that matters to the project beyond the
  current chat (a decision, a constraint, a convention, a source of truth),
  record it in the appropriate project file — this file, or a document under
  `docs/`. Details that only concern the current conversation stay in the chat.
- "Can you ..." is a real question about capability, not a rhetorical request
  to do the thing. Answer whether it is possible and how; do not start doing it
  unless separately asked.

## Parameters are user-controlled and must remain unchanged

- Never change sample, material, device, model, fitting, or numerical parameters.
  This includes notebook settings, parameter files, defaults, and derived
  parameter assignments such as gamma.
- This prohibition also applies to temporary experiments, diagnostic scripts,
  alternate configurations, and in-memory overrides. Do not substitute values
  from publications, tune parameters to improve agreement, or change smoothing
  or grid settings.
- Investigate and correct formulas, data flow, and implementation with the
  existing parameters held fixed. A request to investigate or fix a fit does
  not authorize parameter changes.
- Only an explicit user instruction changing a specific parameter overrides
  this rule for that parameter. Do not repeatedly ask to change parameters.
- Preserve existing user edits. When reverting your work, revert only your own
  changes; do not restore entire files from Git over pre-existing modifications.

## Derive the implementation from basic principles

- Build the calculation from the governing physical equations, with explicit
  assumptions, density definitions, units, signs, and boundary conditions.
  Document the derivation connecting those equations to the implemented steps.
- The user's XLSX workbook is their independent validation tool. Do not read,
  copy, reverse-engineer, or reproduce its formulas or cached results to build
  or repair the script. Do not use exported example results as substitutes for
  deriving the calculation either.
- Use the supplied PDFs as scientific references, checking the original
  equations visually when extracted text is ambiguous. Resolve apparent
  contradictions through derivation and state any remaining uncertainty;
  do not silently choose whichever interpretation makes a curve look right.
- Validate with dimensional consistency, conservation relations, physical
  limits, and independently derived analytical cases. A visually improved fit
  or a test that repeats the implementation does not establish correctness.
- Keep all existing parameters fixed under the rule above. Do not compensate
  for an implementation problem with parameter tuning, arbitrary multipliers,
  clipping, or preprocessing changes.
