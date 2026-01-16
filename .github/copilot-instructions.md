## Engineering and workflow conventions

## Workflow
Use the Makefile : `make build`, `make test`, which internally use the `stack` CLI.

### Clarity before code
Do not jump into coding new features or fixing bugs without first fully understanding the requirements and implications.
Always discuss with the user or team to clarify any ambiguities before starting implementation, and ask follow up questions if needed.

### NO SHORTCUTS
Never leave incomplete functions or stubs in the codebase. Every function must be fully implemented or explicitly marked as unimplemented with a clear error message. 
Never delete or comment out tests to bypass failures; instead, investigate and resolve the underlying issues. 
Never, ever, use "In a real implementation," or similar phrases to justify incomplete code or post-hoc defaults.
Do /not/ mark functions as TODO or FIXME without prior approval from the user. Scope down a feature if needed, and tell the user what's missing.

### Testing and validation
A feature is not complete until it has comprehensive tests that validate its correctness and robustness.
Meaningful tests rather than coverage: Aim for meaningful, high-quality tests that thoroughly validate core functionality and edge cases.
Do not test trivial properties of data structures, but aim for domain-specific validation.


### Debugging and troubleshooting
When investigating or debugging, do not skip or modify tests, or mark them as pending, without receiving prior approval from the user.


### Module Imports
All imports must be explicit. If a module does not export an explicit list, import only the needed symbols. When importing an external library, import only the required functions/types.

### Errors and exceptions
* Never use `error` or `undefined`
* functions with alternative return types should use sum types (e.g. Maybe, Either, or custom ones)
* All error messages should be informative and include context about the failure that can be used to fix the issue.

