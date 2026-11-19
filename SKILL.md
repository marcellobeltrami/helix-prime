---
name: rust-software-development
description: Write, review, refactor, and optimize Rust software with bug prevention, memory optimization, high readability, and thorough comments. Use when the user asks to write Rust code, build Rust crates or services, audit Rust for bugs or memory usage, or clean up and comment existing Rust implementations.
---
# rust-software-development

Produce idiomatic, bug-resistant, memory-conscious, and highly readable Rust code with clear documentation and inline comments.

## When to Use

- Writing new Rust libraries, CLI tools, services, or algorithms from scratch.
- Refactoring or reviewing existing Rust code for correctness, edge cases, or clarity.
- Reducing heap allocations, improving cache locality, or tuning memory layout in Rust.
- Adding docstrings and inline comments to clarify invariants, ownership, or non-obvious logic.

## Steps

1. **Model types and invariants explicitly.** Encode valid states in enums and newtypes rather than loose booleans, sentinel values, or unchecked primitives. Document every struct or enum invariant in a doc comment so callers and maintainers know what must hold true.

2. **Eliminate avoidable bugs and panics.** Prefer `Result` and `Option` over `.unwrap()` or `.expect()` in production paths, and handle edge cases such as empty inputs, integer overflow (`checked_*`, `saturating_*`, or explicit bounds), and boundary indices. Keep `unsafe` blocks minimal, isolated behind safe abstractions, and accompanied by a `// SAFETY:` comment explaining why every soundness requirement is met.

3. **Apply memory optimizations where they add real value.** Minimize unnecessary heap allocations and cloning by borrowing (`&str`, `&[T]`, `Cow<'_, str>`) where ownership is not needed, pre-allocating collections with `Vec::with_capacity` when size is known or bounded, and choosing compact types or `Box<[T]>` over `Vec<T>` when capacity will not grow. Order struct fields from largest to smallest alignment when padding is noticeable, and prefer stack or arena/iterator pipelines over intermediate allocations when benchmarks or hot paths justify it.

4. **Prioritize readability over cleverness.** Write clear, direct control flow instead of dense nested closures, chained iterator tricks, or cryptic index gymnastics. Use descriptive variable and function names, early returns to reduce nesting, and straightforward loops when a loop is easier to read than a complex combinator chain.

5. **Comment intent, invariants, and trade-offs.** Add `///` doc comments on public modules, types, functions, and error variants explaining what they do and any edge cases or failure modes. Add `//` inline comments before non-obvious steps, memory-layout choices, or bug-prevention guards, explaining why the code is written that way rather than restating obvious syntax.

6. **Verify before finishing.** Run `cargo check` and `cargo test` (and `cargo clippy` when available) to confirm compilation, lint cleanliness, and test coverage of normal cases, edge cases, and error paths.

## Gotchas

- Do not sacrifice readability for micro-optimizations unless a hot path or memory constraint benefits from them; clear code that can be audited is easier to keep bug-free.
- Avoid over-using `.clone()` just to satisfy the borrow checker; restructure lifetimes or borrowing scopes first, and comment any intentional clone explaining why ownership is required.
- Do not write comments that merely narrate what a single obvious line already says; reserve comments for invariants, edge cases, non-obvious algorithms, and memory or safety rationale.
- Remember that `Vec::with_capacity(n)` reserves capacity for `n` elements, not length `n`; do not index into an empty pre-allocated vector before pushing or initializing elements.
