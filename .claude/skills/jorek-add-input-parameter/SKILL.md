---
name: jorek-add-input-parameter
description: Wire up a new JOREK namelist input parameter (a new physics flag, switch, or numerical constant read from the input file) end-to-end, following this codebase's established 5-step pattern (declare, namelist, default, log, broadcast). Use this whenever asked to add a new input parameter, add a namelist variable, add a new flag/switch to JOREK, expose a new constant as user input, or make an existing hardcoded value configurable from the input file. Also use it to check/review an existing in-progress addition of such a parameter for completeness (e.g. "did I wire up this parameter correctly?").
---

# Adding a JOREK input parameter

JOREK input parameters (anything read from the `in1` namelist in a user's
input file) are **not** wired up in one place. They touch five files, and
missing any one of them causes a silent bug: the parameter either doesn't
compile, isn't readable from the input file, silently keeps a wrong default,
doesn't show up in the log (making runs hard to reproduce/debug), or works on
MPI rank 0 but not on other ranks (because it was never broadcast).

Go through all five steps below, in order, for every new parameter. Treat
this as a checklist, not a suggestion — half-wiring a parameter is worse than
not adding it, because it fails silently in parallel runs or logs the wrong
value.

A complete worked example lives in git history / current code for the
`freeze_variable` parameter — grep for `freeze_variable` across the five
files named below to see a real instance of every step. It's a derived type
grouping several logical switches plus one real parameter, so it illustrates
both the scalar and derived-type cases.

## The five steps

### 1. Declare it — `models/phys_module.f90`

Add the variable near physically related parameters, with a doxygen `!<`
comment describing its meaning/units:

```fortran
real*8  :: implicit_heat_source !< Choose = 1.d0 to fully switch on the implicit heat source for numerical stabilization
```

**If you're adding a group of related flags** (e.g. several booleans that
belong to one physics feature), don't declare them as loose variables —
define a derived type and a single instance, following the precedent of
`type_freeze_variable`, `type_valve`, `type_bcs`:

```fortran
type :: type_freeze_variable
  logical :: var_psi    !< freeze the psi/mag pot 3 variable
  ...
  real*8  :: zbig       !< large value multiplying the diagonal matrix element of a frozen variable
end type type_freeze_variable

type (type_freeze_variable) :: freeze_variable
```

Use a **scalar** derived type (one global instance) when the feature is a
single on/off switch for the whole run — that's what `type_freeze_variable`
is. Use an **array** of a derived type (e.g. `type_valve, dimension(n_valves_max) :: valves`,
or `type_bcs, dimension(max_bnd_types) :: bcs`) when the setting is
per-boundary, per-group, or otherwise repeated N times.

Don't confuse this file with `models/<model>/mod_model_settings.f90`, which
is a *different* per-model file defining the `var_psi`, `var_u`, ... integer
index parameters (positions into `n_var`-sized arrays like `amat`). Those are
compile-time model configuration, not input parameters — if your new
parameter needs to reference per-variable indices (like `freeze_variable`
does), it reads from `mod_model_settings.f90` but is *declared* in
`phys_module.f90`.

### 2. Add it to the namelist — `models/<model>/initialise_parameters.f90`

Every model directory (`models/model600/`, `models/model307/`, ...) has its
own `initialise_parameters.f90` with a `namelist /in1/ ...` statement listing
every readable parameter. Add your variable there, in whichever model(s) it's
actually relevant to:

```fortran
namelist /in1/  tstep, nstep, ...                                   &
                T_min_neg, rho_min_neg, implicit_heat_source,       &
                freeze_variable,                                    &
                ...
```

A derived type can be listed by its bare name — Fortran namelists support
this, and the user's input file then sets individual fields with
`freeze_variable%var_psi = .true.` syntax.

**Don't assume every model needs it.** Many parameters in `phys_module.f90`
are global/shared but only meaningful for one model's physics (e.g.
`freeze_variable` is declared in the shared module but only appears in
model600's namelist, since the feature only makes sense there). Check which
model(s) actually use the new parameter before deciding where to add the
namelist entry — grep other models for the same parameter name if unsure.

### 3. Set a default — `models/preset_parameters.f90`, `subroutine preset_parameters`

Every parameter needs an explicit default, set before the namelist is read,
so a run that doesn't mention it in the input file still behaves predictably:

```fortran
implicit_heat_source = 0.d0

freeze_variable%var_psi    = .false.
...
freeze_variable%zbig       = 1.d12
```

Logicals almost always default to `.false.` (feature off) unless there's a
specific physics reason otherwise — check for a comment or ask if unsure.
Reals need a physically meaningful default, not just `0.d0` by reflex (e.g.
`zbig = 1.d12` because it has to dominate the diagonal matrix element it's
multiplying).

### 4. Log it — `models/mod_log_params.f90`, `subroutine log_parameters`

So that every run's stdout log records what was actually used (essential for
reproducing/debugging a run later). Use the existing format constants defined
near the top of the file:

- `REAL_FMT` for `real*8`
- `LOGI_FMT` for `logical`
- `CHAR_FMT` for `character`

```fortran
write(*,REAL_FMT) 'implicit_heat_source  ', implicit_heat_source
write(*,LOGI_FMT) 'freeze_variable%var_psi   ', freeze_variable%var_psi
...
```

Place it near the related parameters so the log output stays grouped and
readable.

### 5. Broadcast it — `communication/broadcast_phys.f90`

Parameters are only read on MPI rank 0 (from the input file), so they must be
packed into a buffer and broadcast to every other rank, or those ranks will
silently keep the preset default forever. This file has two halves in the
same subroutine: a pack block (`if (my_id .eq. 0) then`) and an unpack block
(`if (my_id .ne. 0) then`) further down. Add one `MPI_PACK` call and one
matching `MPI_UNPACK` call, **in the same relative order** in both blocks —
order must match exactly since unpacking reads the buffer sequentially.

```fortran
! pack side
call MPI_PACK(implicit_heat_source,   1,MPI_REAL8,buffer,bufsize,position,MPI_COMM_WORLD,ierr)
call MPI_PACK(freeze_variable%var_psi,1,MPI_LOGICAL,buffer,bufsize,position,MPI_COMM_WORLD,ierr)

! unpack side, same order
call MPI_UNPACK(buffer,bufsize,position,implicit_heat_source,   1,MPI_REAL8,MPI_COMM_WORLD,ierr)
call MPI_UNPACK(buffer,bufsize,position,freeze_variable%var_psi,1,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
```

Type mapping: `MPI_LOGICAL` for `logical`, `MPI_REAL8` for `real*8`,
`MPI_INTEGER` for `integer`, `MPI_CHARACTER` for `character` (with the string
length as the count argument, e.g. `n_valves_max*4` for an array of
4-character strings).

JOREK's convention here is always **field-by-field** packing of derived
types — never pack a whole derived type in one `MPI_PACK` call, even though
Fortran would technically allow it for a sequence-associated type. This
keeps every broadcast's wire format explicit and greppable.

`bufsize` is a fixed, generous constant (currently `2500000`) — you don't
need to recompute or account for buffer size when adding a parameter.

## After wiring it up

Grep for the parameter name across all five files to confirm every step
landed:

```
grep -rn "<param_name>" models/phys_module.f90 models/preset_parameters.f90 \
  models/mod_log_params.f90 communication/broadcast_phys.f90 \
  models/<model>/initialise_parameters.f90
```

Compiling Fortran on this shared login node needs the user's explicit
go-ahead — don't compile to "double check" without asking first.
