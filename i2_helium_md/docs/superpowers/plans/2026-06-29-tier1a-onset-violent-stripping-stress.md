# Tier-1a Onset-Violent Stripping Stress Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a separate Tier-1a stress-test schedule family that strips the He shell violently at `t=0.5 ps` to endpoints `n={14,2,1,0}`, using continuous velocity and separate stress reporting.

**Architecture:** Reuse the existing Tier-1a schedule/driver/reporting seams. Add a schedule-compatible one-event `onset_strip` builder, generalize the continuous-velocity shed primitive and `shed_step` from one removed He to `n_before - n_after`, then add stress-specific run naming, generation, scoring, and plotting scripts that stay separate from the physical Tier-1a RMSE table.

**Tech Stack:** Python 3.14, NumPy, dataclasses, pytest, existing `i2_helium_md` simulation/postprocess modules, existing Tier-1a script patterns.

---

## File Structure

- Modify `i2_helium_md/physics/shell_schedule.py`: add `build_onset_strip_schedule(...)`; keep `build_shell_schedule(...)` unchanged.
- Modify `i2_helium_md/physics/mass_jump.py`: add `n_removed` support to continuous-velocity shedding; keep cold-shed bound helpers unchanged.
- Modify `i2_helium_md/config.py`: add explicit `anchor_n_final` endpoint metadata for stress schedules.
- Modify `i2_helium_md/simulation/ion_propagation_step.py`: make `shed_step(...)` compute `n_removed` from each event.
- Modify `scripts/tier1a_common.py`: add stress constants, config builder, and stress run-tag helpers.
- Create `scripts/gen_tier1a_stress_runs.py`: generate fixed null plus four onset-strip stress runs.
- Create `scripts/post_processing/tier1a_stress_table.py`: score stress runs separately and plot fixed plus one selected stress case.
- Modify `tests/test_shell_schedule.py`: add onset-strip schedule tests.
- Modify `tests/test_mass_jump.py`: add batch continuous-velocity shed tests.
- Modify `tests/test_ion_variable_mass.py`: add batch-event `shed_step` test.
- Modify `tests/test_ion_drag_smoke.py`: add run-level stress smoke test.
- Modify `tests/test_tier1a_scripts.py`: add stress naming/scorer/plot tests.
- Modify `TIER1A_IMPLEMENTATION_PLAN.md`: add completed diagnostic stress slice when implementation is complete.
- Modify `drag_migration_log_tier1a.md`: add delivery record when implementation is complete.

---

### Task 1: Add Onset-Strip Schedule Tests

**Files:**
- Modify: `tests/test_shell_schedule.py`
- Later implementation target: `i2_helium_md/physics/shell_schedule.py`

- [ ] **Step 1: Write failing tests for one-event onset-strip schedule**

Add imports if missing:

```python
import pytest

from i2_helium_md.physics.shell_schedule import (
    ANCHOR_N_START,
    build_onset_strip_schedule,
    complex_mass_amu,
)
```

Add tests:

```python
class TestOnsetStripSchedule:
    @pytest.mark.parametrize("n_final,n_removed", [(14, 7), (2, 19), (1, 20), (0, 21)])
    def test_one_event_at_onset_with_requested_endpoint(self, n_final, n_removed):
        sched = build_onset_strip_schedule(t_strip_ps=0.5, n_final=n_final)

        assert sched.t_star_ps == pytest.approx(0.5)
        assert sched.crossing_fraction == pytest.approx(1.0)
        assert len(sched.events) == 1

        event = sched.events[0]
        assert event.index == 1
        assert event.time_ps == pytest.approx(0.5)
        assert event.n_before == ANCHOR_N_START
        assert event.n_after == n_final
        assert event.n_before - event.n_after == n_removed
        assert event.mass_before_amu == pytest.approx(complex_mass_amu(21))
        assert event.mass_after_amu == pytest.approx(complex_mass_amu(n_final))

    @pytest.mark.parametrize("n_final", [14, 2, 1, 0])
    def test_n_of_t_jumps_directly_at_onset(self, n_final):
        sched = build_onset_strip_schedule(t_strip_ps=0.5, n_final=n_final)

        assert sched.n_of_t(0.499999) == 21
        assert sched.n_of_t(0.5) == n_final
        assert sched.n_of_t(30.0) == n_final

    def test_n_bar_matches_step_for_stress_schedule(self):
        sched = build_onset_strip_schedule(t_strip_ps=0.5, n_final=2)

        assert sched.n_bar(0.0) == pytest.approx(21.0)
        assert sched.n_bar(0.5) == pytest.approx(2.0)
        assert sched.n_bar(30.0) == pytest.approx(2.0)

    @pytest.mark.parametrize("bad_n_final", [-1, 21, 22])
    def test_rejects_invalid_endpoint(self, bad_n_final):
        with pytest.raises(ValueError, match="n_final"):
            build_onset_strip_schedule(t_strip_ps=0.5, n_final=bad_n_final)

    def test_rejects_negative_strip_time(self):
        with pytest.raises(ValueError, match="t_strip_ps"):
            build_onset_strip_schedule(t_strip_ps=-0.1, n_final=0)
```

- [ ] **Step 2: Run schedule tests and verify failure**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_shell_schedule.py -q
```

Expected: fail during import with `cannot import name 'build_onset_strip_schedule'`.

- [ ] **Step 3: Implement `build_onset_strip_schedule`**

In `i2_helium_md/physics/shell_schedule.py`, add this class method override to `ShellSchedule`:

```python
    def n_bar(self, t_ps) -> np.ndarray | float:
        """Continuous or stress shell-count guide [dimensionless]."""
        if len(self.events) == 1 and self.crossing_fraction == 1.0:
            t = np.asarray(t_ps, dtype=float)
            event = self.events[0]
            n = np.where(t < event.time_ps, float(event.n_before), float(event.n_after))
            return float(n) if np.ndim(t_ps) == 0 else n
        t = np.asarray(t_ps, dtype=float)
        ts = self.t_star_ps
        slope1 = (ANCHOR_N_START - ANCHOR_N_MID) / (ANCHOR_T_MID_PS - ts)
        slope2 = (ANCHOR_N_MID - ANCHOR_N_END) / (ANCHOR_T_END_PS - ANCHOR_T_MID_PS)
        n = np.select(
            [t <= ts, t <= ANCHOR_T_MID_PS, t <= ANCHOR_T_END_PS],
            [
                float(ANCHOR_N_START),
                ANCHOR_N_START - slope1 * (t - ts),
                ANCHOR_N_MID - slope2 * (t - ANCHOR_T_MID_PS),
            ],
            default=float(ANCHOR_N_END),
        )
        return float(n) if np.ndim(t_ps) == 0 else n
```

Then add the builder below `build_shell_schedule(...)`:

```python
def build_onset_strip_schedule(
    *,
    t_strip_ps: float = 0.5,
    n_final: int,
) -> ShellSchedule:
    """Build a one-event onset-violent stripping stress schedule.

    This is a diagnostic stress schedule, not a TDDFT-anchored physical loss curve.
    It strips directly from n=21 to ``n_final`` at ``t_strip_ps``.
    """
    if not np.isfinite(t_strip_ps) or t_strip_ps < 0.0:
        raise ValueError(f"t_strip_ps must be finite and >= 0; got {t_strip_ps!r}.")
    if int(n_final) != n_final or not (0 <= int(n_final) < ANCHOR_N_START):
        raise ValueError(
            f"n_final must be an integer in [0, {ANCHOR_N_START - 1}], got {n_final!r}."
        )

    n_after = int(n_final)
    mass_before = complex_mass_amu(ANCHOR_N_START)
    mass_after = complex_mass_amu(n_after)
    event = ShedEvent(
        index=1,
        crossing=float(n_after),
        time_ps=float(t_strip_ps),
        n_before=ANCHOR_N_START,
        n_after=n_after,
        mass_before_amu=mass_before,
        mass_after_amu=mass_after,
        kick_factor=mass_before / mass_after,
    )
    return ShellSchedule(
        t_star_ps=float(t_strip_ps),
        crossing_fraction=1.0,
        events=(event,),
    )
```

- [ ] **Step 4: Run schedule tests and verify pass**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_shell_schedule.py -q
```

Expected: all `tests/test_shell_schedule.py` tests pass.

- [ ] **Step 5: Commit Task 1**

```powershell
git add i2_helium_md/physics/shell_schedule.py tests/test_shell_schedule.py
git commit -m "Add onset strip stress schedule"
```

---

### Task 2: Generalize Continuous-Velocity Shed To Batch Removal

**Files:**
- Modify: `tests/test_mass_jump.py`
- Modify: `i2_helium_md/physics/mass_jump.py`

- [ ] **Step 1: Write failing batch shed tests**

In `tests/test_mass_jump.py`, add:

```python
class TestBatchContinuousVelocityShed:
    @pytest.mark.parametrize("n_after", [14, 2, 1, 0])
    def test_batch_velocity_unchanged_mass_drop_and_energy(self, n_after):
        m = complex_mass_amu(21)
        n_removed = 21 - n_after
        res = continuous_velocity_shed(V_MINUS, m, n_removed=n_removed)

        np.testing.assert_array_equal(res.v_plus, V_MINUS)
        assert res.m_plus_amu == pytest.approx(complex_mass_amu(n_after), abs=1e-12)
        speed_sq = float(V_MINUS @ V_MINUS)
        assert res.dE_mass_transfer == pytest.approx(
            0.5 * n_removed * MASS_HE_AMU * speed_sq,
            rel=1e-12,
        )

    def test_batch_total_momentum_and_ke_with_removed_co_moving_he_conserved(self):
        m = complex_mass_amu(21)
        n_removed = 21
        res = continuous_velocity_shed(V_MINUS, m, n_removed=n_removed)
        removed_mass = n_removed * MASS_HE_AMU

        momentum_before = m * V_MINUS
        momentum_after = res.m_plus_amu * res.v_plus + removed_mass * V_MINUS
        np.testing.assert_allclose(momentum_after, momentum_before, rtol=0.0, atol=1e-12)

        speed_sq = float(V_MINUS @ V_MINUS)
        ke_before = 0.5 * m * speed_sq
        ke_after = 0.5 * res.m_plus_amu * speed_sq + 0.5 * removed_mass * speed_sq
        assert ke_after == pytest.approx(ke_before, rel=1e-12)

    @pytest.mark.parametrize("bad_n_removed", [0, -1, 100])
    def test_rejects_invalid_batch_removed_count(self, bad_n_removed):
        with pytest.raises(ValueError, match="n_removed"):
            continuous_velocity_shed(V_MINUS, complex_mass_amu(21), n_removed=bad_n_removed)
```

In `tests/test_ion_variable_mass.py`, extend `TestVectorizedContinuousVelocityShed`:

```python
    def test_batch_vectorized_energy_uses_removed_count(self):
        m = complex_mass_amu(21)
        vx = np.array([2.0])
        vy = np.array([0.0])
        vz = np.array([0.0])
        vxp, vyp, vzp, m_plus, dE = continuous_velocity_shed_components(
            vx, vy, vz, m, n_removed=21,
        )

        np.testing.assert_array_equal(vxp, vx)
        np.testing.assert_array_equal(vyp, vy)
        np.testing.assert_array_equal(vzp, vz)
        assert m_plus == pytest.approx(complex_mass_amu(0))
        assert dE[0] == pytest.approx(0.5 * 21 * MASS_HE_AMU * 4.0)
```

Add `MASS_HE_AMU` import to `tests/test_ion_variable_mass.py`:

```python
from i2_helium_md.physics.constants import MASS_HE_AMU, U
```

- [ ] **Step 2: Run targeted tests and verify failure**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_mass_jump.py tests/test_ion_variable_mass.py -q
```

Expected: fail with `unexpected keyword argument 'n_removed'`.

- [ ] **Step 3: Implement batch support in scalar and vector helpers**

In `i2_helium_md/physics/mass_jump.py`, add the mass-agnostic count validator below.
It must not infer shell count from iodine mass; shell-count legality is owned by the
schedule layer. This helper validates only the count type, positivity, and positive
post-shed mass.

```python
def _check_removed_count(n_removed: int, m_minus_amu: float, m_he_amu: float) -> int:
    """Return validated integer removed-He count."""
    if isinstance(n_removed, (bool, np.bool_)) or not isinstance(n_removed, (int, np.integer)):
        raise ValueError(f"n_removed must be an integer; got {n_removed!r}.")
    n = int(n_removed)
    if n <= 0:
        raise ValueError(f"n_removed must be > 0; got {n_removed!r}.")
    if not (m_minus_amu - n * m_he_amu > 0.0):
        raise ValueError(
            f"n_removed={n} removes too much mass from m_minus_amu={m_minus_amu!r}."
        )
    return n
```

Change `continuous_velocity_shed(...)` signature:

```python
def continuous_velocity_shed(
    v_minus,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
) -> ShedResult:
```

Replace its mass/energy body with:

```python
    n = _check_removed_count(n_removed, m_minus_amu, m_he_amu)
    m_plus = m_minus_amu - n * m_he_amu
    speed_sq = float(np.sum(v ** 2))
    dE_mass_transfer = 0.5 * n * m_he_amu * speed_sq
```

Change `continuous_velocity_shed_components(...)` signature:

```python
def continuous_velocity_shed_components(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
```

Replace its mass/energy body with:

```python
    n = _check_removed_count(n_removed, m_minus_amu, m_he_amu)
    m_plus = m_minus_amu - n * m_he_amu
    speed_sq = vxf ** 2 + vyf ** 2 + vzf ** 2
    dE_mass_transfer = 0.5 * n * m_he_amu * speed_sq
```

Change `apply_shed(...)` signature:

```python
def apply_shed(
    v_minus,
    m_minus_amu: float,
    *,
    mode: ShedMode,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
) -> ShedResult:
```

Change anchored delegation:

```python
        return continuous_velocity_shed(
            v_minus,
            m_minus_amu,
            m_he_amu=m_he_amu,
            n_removed=n_removed,
        )
```

- [ ] **Step 4: Run targeted tests and verify pass**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_mass_jump.py tests/test_ion_variable_mass.py -q
```

Expected: both files pass.

- [ ] **Step 5: Commit Task 2**

```powershell
git add i2_helium_md/physics/mass_jump.py tests/test_mass_jump.py tests/test_ion_variable_mass.py
git commit -m "Support batch continuous velocity shedding"
```

---

### Task 3: Wire Batch Events Through `shed_step`

**Files:**
- Modify: `tests/test_ion_variable_mass.py`
- Modify: `tests/test_ion_drag_smoke.py`
- Modify: `i2_helium_md/simulation/ion_propagation_step.py`
- Modify: `i2_helium_md/simulation/ion.py`

- [ ] **Step 1: Write failing driver unit test for batch event**

In `tests/test_ion_variable_mass.py`, import the stress schedule:

```python
from i2_helium_md.physics.shell_schedule import (
    build_onset_strip_schedule,
    build_shell_schedule,
    complex_mass_amu,
)
```

Add to `TestShedStep`:

```python
    def test_batch_event_updates_mass_and_keeps_velocity(self):
        sched = build_onset_strip_schedule(t_strip_ps=0.5, n_final=0)
        st = _state(0.495, vx=3.0, vy=4.0, mass_amu=complex_mass_amu(21))
        new, idx = shed_step(st, sched, 0, dt=0.01)

        assert idx == 1
        assert new.mass_kg[0] == pytest.approx(complex_mass_amu(0) * U)
        assert new.vx[0] == pytest.approx(3.0)
        assert new.vy[0] == pytest.approx(4.0)
        assert new.vz[0] == pytest.approx(0.0)
        assert np.all(new.E_mass_transfer_eV > 0.0)
```

- [ ] **Step 2: Run unit tests and verify failure**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_ion_variable_mass.py -q
```

Expected: batch test fails because `shed_step` removes only one He.

- [ ] **Step 3: Modify `shed_step` to pass removed count**

In `i2_helium_md/simulation/ion_propagation_step.py`, replace the call with:

```python
    n_removed = event.n_before - event.n_after
    vx_p, vy_p, vz_p, m_plus_amu, dE_amu = continuous_velocity_shed_components(
        state.vx,
        state.vy,
        state.vz,
        event.mass_before_amu,
        m_he_amu=m_he_amu,
        n_removed=n_removed,
    )
```

Update the `shed_step` docstring sentence:

```python
    velocities unchanged, drop the scheduled number of He atoms from the uniform
    complex mass, book the per-atom co-moving-He kinetic energy into
    ``E_mass_transfer_eV``, and advance the pointer.
```

- [ ] **Step 4: Add config endpoint field and driver selection for stress schedules**

In `i2_helium_md/config.py`, extend `AnchorMode`:

```python
AnchorMode = Literal["time", "onset_strip"]
```

Add the explicit stress endpoint field next to `anchor_mode`:

```python
    anchor_mode: AnchorMode = "time"      # time-anchored; "onset_strip" is diagnostic stress
    anchor_n_final: int = 14              # endpoint for onset_strip stress only
```

In `SimConfig.validate(...)`, add this guard after the anchor-mode validation:

```python
        if self.anchor_mode == "onset_strip":
            if int(self.anchor_n_final) != self.anchor_n_final:
                raise ValueError(
                    f"anchor_n_final must be an integer; got {self.anchor_n_final!r}."
                )
            if not (0 <= int(self.anchor_n_final) < 21):
                raise ValueError(
                    f"anchor_n_final must be in [0, 20] for onset_strip; "
                    f"got {self.anchor_n_final!r}."
                )
```

In `i2_helium_md/simulation/ion.py`, find the schedule construction. Replace the anchored-only call with:

```python
        if cfg.anchor_mode == "time":
            schedule = build_shell_schedule(cfg.t_star_ps)
        elif cfg.anchor_mode == "onset_strip":
            schedule = build_onset_strip_schedule(
                t_strip_ps=cfg.t_star_ps,
                n_final=cfg.anchor_n_final,
            )
        else:
            raise NotImplementedError(f"unsupported anchor_mode={cfg.anchor_mode!r}")
```

At the top of `ion.py`, extend the existing shell-schedule import:

```python
from ..physics.shell_schedule import build_onset_strip_schedule, build_shell_schedule
```

If `build_shell_schedule` is already imported, extend the existing import rather than duplicating it.

- [ ] **Step 5: Run driver unit tests**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_ion_variable_mass.py tests/test_ion_drag_smoke.py -q
```

Expected: existing anchored tests still pass; batch unit test passes.

- [ ] **Step 6: Commit Task 3**

```powershell
git add i2_helium_md/config.py i2_helium_md/simulation/ion.py i2_helium_md/simulation/ion_propagation_step.py tests/test_ion_variable_mass.py tests/test_ion_drag_smoke.py
git commit -m "Wire onset strip stress schedules through ion driver"
```

---

### Task 4: Add Stress Config And Naming Helpers

**Files:**
- Modify: `scripts/tier1a_common.py`
- Modify: `tests/test_tier1a_scripts.py`

- [ ] **Step 1: Write failing helper tests**

In `tests/test_tier1a_scripts.py`, add:

```python
def test_build_onset_strip_cfg_sets_stress_fields():
    from scripts.tier1a_common import build_onset_strip_cfg

    cfg = build_onset_strip_cfg(
        "9A",
        "shared_pure_cubic",
        n_final=2,
        t_strip_ps=0.5,
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
    )

    assert cfg.mass_scenario == "anchored_discrete"
    assert cfg.anchor_mode == "onset_strip"
    assert cfg.anchor_n_final == 2
    assert cfg.t_star_ps == 0.5
    assert cfg.mass_initial_amu == pytest.approx(complex_mass_amu(21))
    assert cfg.allow_inconsistent_mass_pairing is True


def test_tier1a_stress_run_names_are_distinct():
    from scripts.tier1a_common import tier1a_stress_run_dir_name, tier1a_stress_run_tag

    assert (
        tier1a_stress_run_tag(n_final=0, t_strip_ps=0.5)
        == "tier1a_stress_onset_strip_n0_t0.5"
    )
    assert (
        tier1a_stress_run_dir_name("9A", "shared_pure_cubic", 50, n_final=0, t_strip_ps=0.5)
        == "9A_drag_shared_pure_cubic_N50_tier1a_stress_onset_strip_n0_t0.5"
    )
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier1a_scripts.py -q
```

Expected: import failure for `build_onset_strip_cfg`.

- [ ] **Step 3: Implement stress helpers**

In `scripts/tier1a_common.py`, add constants:

```python
TIER1A_STRESS_T_STRIP_PS = 0.5
TIER1A_STRESS_N_FINAL_VALUES: tuple[int, ...] = (14, 2, 1, 0)
```

Add builder:

```python
def build_onset_strip_cfg(
    case: str,
    variant: str,
    *,
    n_final: int,
    t_strip_ps: float,
    num_molecules: int,
    ion_time_ps: float,
    dt_ion_ps: float,
    seed: int,
    coeff_overrides: Optional[Mapping[str, float]] = None,
    e_bind_override: Optional[float] = None,
) -> SimConfig:
    """Build a diagnostic Tier-1a onset-strip stress config."""
    fixed_cfg = build_drag_cfg(
        case,
        variant,
        num_molecules=num_molecules,
        ion_time_ps=ion_time_ps,
        dt_ion_ps=dt_ion_ps,
        seed=seed,
        coeff_overrides=coeff_overrides,
        e_bind_override=e_bind_override,
    )
    cfg = replace(
        fixed_cfg,
        mass_scenario="anchored_discrete",
        t_star_ps=float(t_strip_ps),
        anchor_mode="onset_strip",
        anchor_n_final=int(n_final),
        coulomb_available_eV=0.80,
        allow_inconsistent_mass_pairing=True,
        mass_initial_amu=complex_mass_amu(21),
    )
    cfg.validate()
    return cfg
```

Add naming:

```python
def tier1a_stress_run_tag(*, n_final: int, t_strip_ps: float) -> str:
    """Return the run-tag suffix for onset-strip stress runs."""
    return f"tier1a_stress_onset_strip_n{int(n_final)}_t{float(t_strip_ps):.1f}"


def tier1a_stress_run_dir_name(
    case: str,
    variant: str,
    n: int,
    *,
    n_final: int,
    t_strip_ps: float,
) -> str:
    """Return the Tier-0-style run directory basename for an onset-strip stress run."""
    return run_dir_name(
        case,
        variant,
        n,
        run_tag=tier1a_stress_run_tag(n_final=n_final, t_strip_ps=t_strip_ps),
    )
```

Update `__all__` with the new constants and functions.

- [ ] **Step 4: Run helper tests and verify pass**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier1a_scripts.py -q
```

Expected: pass with the existing intentional pairing warning.

- [ ] **Step 5: Commit Task 4**

```powershell
git add scripts/tier1a_common.py tests/test_tier1a_scripts.py
git commit -m "Add Tier-1a onset strip stress helpers"
```

---

### Task 5: Add Stress Generator Script

**Files:**
- Create: `scripts/gen_tier1a_stress_runs.py`

- [ ] **Step 1: Create generator script**

Create `scripts/gen_tier1a_stress_runs.py`:

```python
"""Generate Tier-1a onset-violent stripping stress runs.

Produces one fixed null and four diagnostic onset-strip stress run directories
under ``data/runs``. These are sensitivity tests, not TDDFT-anchored physical
Tier-1a runs.
"""

from __future__ import annotations

from pathlib import Path
import sys


CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 50

ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604
T_STRIP_PS = 0.5
N_FINAL_VALUES = [14, 2, 1, 0]

COEFF_OVERRIDES: dict[str, float] = {}
E_BIND_OVERRIDE = None


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier0_common import build_drag_cfg  # noqa: E402
from scripts.tier1a_common import (  # noqa: E402
    build_onset_strip_cfg,
    tier1a_run_dir_name,
    tier1a_stress_run_dir_name,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one run directory using the standard pipeline."""
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(
        f"[{label}] form={cfg.drag_form} "
        f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
        f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV "
        f"scenario={cfg.mass_scenario} anchor_mode={cfg.anchor_mode} "
        f"t_strip={cfg.t_star_ps:.3f} ps"
    )
    print(f"[{label}] neutral propagation ...")
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

    print(f"[{label}] ion propagation ...")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

    print(
        f"[{label}] done -> {run_dir} "
        f"(neutral {neutral.time_ps.size} steps, ion {ion.time_ps.size} steps)"
    )


def main() -> int:
    fixed_cfg = build_drag_cfg(
        CASE,
        VARIANT,
        num_molecules=N,
        ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS,
        seed=SEED,
        coeff_overrides=COEFF_OVERRIDES,
        e_bind_override=E_BIND_OVERRIDE,
    )
    fixed_dir = PROJECT_ROOT / "data" / "runs" / tier1a_run_dir_name(
        CASE, VARIANT, N, "fixed", None
    )
    _run_one(f"{CASE} {VARIANT} N={N} fixed", fixed_cfg, fixed_dir)

    for n_final in N_FINAL_VALUES:
        cfg = build_onset_strip_cfg(
            CASE,
            VARIANT,
            n_final=n_final,
            t_strip_ps=T_STRIP_PS,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            coeff_overrides=COEFF_OVERRIDES,
            e_bind_override=E_BIND_OVERRIDE,
        )
        run_dir = PROJECT_ROOT / "data" / "runs" / tier1a_stress_run_dir_name(
            CASE, VARIANT, N, n_final=n_final, t_strip_ps=T_STRIP_PS
        )
        _run_one(
            f"{CASE} {VARIANT} N={N} onset-strip n_final={n_final}",
            cfg,
            run_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Smoke import the generator**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m py_compile scripts/gen_tier1a_stress_runs.py
```

Expected: no output, exit code 0.

- [ ] **Step 3: Commit Task 5**

```powershell
git add scripts/gen_tier1a_stress_runs.py
git commit -m "Add Tier-1a stress run generator"
```

---

### Task 6: Add Stress Scorer And Plotter

**Files:**
- Create: `scripts/post_processing/tier1a_stress_table.py`
- Modify: `tests/test_tier1a_scripts.py`

- [ ] **Step 1: Write failing scorer tests**

In `tests/test_tier1a_scripts.py`, add:

```python
def test_score_tier1a_stress_run_emits_stress_columns(tmp_path):
    from scripts.post_processing.tier1a_stress_table import (
        TIER1A_STRESS_TABLE_COLUMNS,
        score_tier1a_stress_run,
    )

    run = RunDirectory(tmp_path / "stress")
    run.save_ion(_tiny_ion_checkpoint())

    row = score_tier1a_stress_run(
        run,
        case="9A",
        variant="shared_pure_cubic",
        n=1,
        run_tag="tier1a_stress_onset_strip_n0_t0.5",
        n_final=0,
        t_strip_ps=0.5,
        hedft=_hedft_reference(),
        smoothed=_smoothed_reference(),
        window_start_ps=1.0,
    )

    assert list(row) == TIER1A_STRESS_TABLE_COLUMNS
    assert row["stress_family"] == "onset_strip"
    assert row["n_final_requested"] == 0
    assert row["t_strip_ps"] == 0.5
    assert row["n_removed"] == 7
```

Add plot test:

```python
def test_build_tier1a_stress_figure_plots_fixed_plus_selected_case():
    import matplotlib

    matplotlib.use("Agg")

    from scripts.post_processing.tier1a_stress_table import (
        Tier1aStressRunRecord,
        build_stress_trajectory_figure,
    )

    ion = _tiny_ion_checkpoint()
    records = [
        Tier1aStressRunRecord(
            label="fixed",
            run_tag="tier1a_fixed",
            n_final=None,
            t_strip_ps=None,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aStressRunRecord(
            label="onset strip n=0",
            run_tag="tier1a_stress_onset_strip_n0_t0.5",
            n_final=0,
            t_strip_ps=0.5,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aStressRunRecord(
            label="onset strip n=2",
            run_tag="tier1a_stress_onset_strip_n2_t0.5",
            n_final=2,
            t_strip_ps=0.5,
            ion=ion,
            cfg=None,
            row={},
        ),
    ]

    fig = build_stress_trajectory_figure(
        records,
        _hedft_reference(),
        _smoothed_reference(),
        window=(1.0, 3.0),
        selected_n_final=0,
        title="Tier-1a stress test",
    )

    labels = {line.get_label() for line in fig.axes[0].lines}
    assert "fixed MD mean |v2|" in labels
    assert "onset strip n=0 MD mean |v2|" in labels
    assert "onset strip n=2 MD mean |v2|" not in labels
```

- [ ] **Step 2: Run script tests and verify failure**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier1a_scripts.py -q
```

Expected: import failure for `scripts.post_processing.tier1a_stress_table`.

- [ ] **Step 3: Create stress table module**

Create `scripts/post_processing/tier1a_stress_table.py` by adapting the Tier-1a table shape:

```python
"""Emit the Tier-1a onset-strip stress table for finished run directories."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 50
T_STRIP_PS = 0.5
N_FINAL_VALUES = [14, 2, 1, 0]
WINDOW_START_PS = 2.67
PLOT_N_FINAL = 0

SAVE_CSV_PATH = None
EXPORT_MEAN_SERIES_DIR = PROJECT_ROOT / "data" / "reference" / "drag" / CASE / "tier1a_stress"
SHOW_FIGURE = True


from scripts.tier1a_common import (  # noqa: E402
    tier1a_run_dir_name,
    tier1a_run_tag,
    tier1a_stress_run_dir_name,
    tier1a_stress_run_tag,
)
from scripts.post_processing.tier0_drag_comparison import (  # noqa: E402
    ensemble_mean_series,
    export_mean_series,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    SmoothedSpeedReference,
    compare_distance,
    compare_speed_to_reference,
    ion_ledger_closure,
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


TIER1A_STRESS_TABLE_COLUMNS = [
    "case",
    "variant",
    "N",
    "run_tag",
    "stress_family",
    "n_final_requested",
    "t_strip_ps",
    "window_start_ps",
    "window_end_ps",
    "n_scored_R",
    "n_scored_v2_smoothed",
    "RMSE_R_raw_A",
    "RMSE_v2_smoothed_Aps",
    "v2_smoothed_mean_ratio",
    "ledger_max_resid_eV",
    "n_removed",
    "n_shell_start",
    "n_shell_end",
]


@dataclass(frozen=True)
class Tier1aStressRunRecord:
    label: str
    run_tag: str
    n_final: int | None
    t_strip_ps: float | None
    ion: IonCheckpoint
    cfg: SimConfig | None
    row: dict[str, Any]


def _load_references(project_root: Path, case: str) -> tuple[HedftTrajectory, SmoothedSpeedReference]:
    hedft = load_hedft_trajectory(project_root / "data" / "reference" / f"{case}_All_Data.csv")
    smoothed = load_smoothed_speed_reference(
        project_root / "data" / "reference" / "drag" / case / "velocity_smoothed" / "cleaned_data_long.csv"
    )
    return hedft, smoothed


def _shell_summary(n_shell: np.ndarray) -> tuple[int, int, int]:
    n_start = int(round(float(n_shell[0, 0])))
    n_end = int(round(float(n_shell[0, -1])))
    return n_start, n_end, max(0, n_start - n_end)


def score_tier1a_stress_run(
    run: RunDirectory | str | Path,
    *,
    case: str,
    variant: str,
    n: int,
    run_tag: str,
    n_final: int,
    t_strip_ps: float,
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window_start_ps: float = WINDOW_START_PS,
) -> dict[str, Any]:
    run_dir = run if isinstance(run, RunDirectory) else RunDirectory(run)
    ion = run_dir.load_ion()
    window_end_ps = float(smoothed.time_ps[-1])
    window = (float(window_start_ps), window_end_ps)
    dist = compare_distance(ion, hedft, window=window)
    v2_smoothed = compare_speed_to_reference(
        ion,
        atom="I2",
        t_ref_ps=smoothed.time_ps,
        ref_speed_Aps=smoothed.speed_Aps,
        window=window,
    )
    closure = ion_ledger_closure(ion)
    n_shell_start, n_shell_end, n_removed = _shell_summary(ion.n_shell)
    return {
        "case": case,
        "variant": variant,
        "N": int(n),
        "run_tag": run_tag,
        "stress_family": "onset_strip",
        "n_final_requested": int(n_final),
        "t_strip_ps": float(t_strip_ps),
        "window_start_ps": window[0],
        "window_end_ps": window[1],
        "n_scored_R": dist.num_overlap_points,
        "n_scored_v2_smoothed": v2_smoothed.num_overlap_points,
        "RMSE_R_raw_A": dist.rmse,
        "RMSE_v2_smoothed_Aps": v2_smoothed.rmse,
        "v2_smoothed_mean_ratio": v2_smoothed.mean_ratio,
        "ledger_max_resid_eV": closure.max_abs_residual_eV,
        "n_removed": n_removed,
        "n_shell_start": n_shell_start,
        "n_shell_end": n_shell_end,
    }


def collect_tier1a_stress_records(
    *,
    project_root: Path = PROJECT_ROOT,
    case: str = CASE,
    variant: str = VARIANT,
    n: int = N,
    n_final_values: Iterable[int] = tuple(N_FINAL_VALUES),
    t_strip_ps: float = T_STRIP_PS,
    window_start_ps: float = WINDOW_START_PS,
) -> tuple[list[Tier1aStressRunRecord], HedftTrajectory, SmoothedSpeedReference]:
    hedft, smoothed = _load_references(project_root, case)
    run_root = project_root / "data" / "runs"
    records: list[Tier1aStressRunRecord] = []

    fixed_tag = tier1a_run_tag("fixed", None)
    fixed_run = RunDirectory(run_root / tier1a_run_dir_name(case, variant, n, "fixed", None))
    fixed_ion = fixed_run.load_ion()
    records.append(
        Tier1aStressRunRecord(
            label="fixed",
            run_tag=fixed_tag,
            n_final=None,
            t_strip_ps=None,
            ion=fixed_ion,
            cfg=fixed_run.load_cfg() if fixed_run.has_cfg() else None,
            row={},
        )
    )

    for n_final in n_final_values:
        run_tag = tier1a_stress_run_tag(n_final=n_final, t_strip_ps=t_strip_ps)
        run = RunDirectory(
            run_root / tier1a_stress_run_dir_name(
                case, variant, n, n_final=n_final, t_strip_ps=t_strip_ps
            )
        )
        ion = run.load_ion()
        row = score_tier1a_stress_run(
            run,
            case=case,
            variant=variant,
            n=n,
            run_tag=run_tag,
            n_final=n_final,
            t_strip_ps=t_strip_ps,
            hedft=hedft,
            smoothed=smoothed,
            window_start_ps=window_start_ps,
        )
        records.append(
            Tier1aStressRunRecord(
                label=f"onset strip n={n_final}",
                run_tag=run_tag,
                n_final=n_final,
                t_strip_ps=t_strip_ps,
                ion=ion,
                cfg=run.load_cfg() if run.has_cfg() else None,
                row=row,
            )
        )
    return records, hedft, smoothed


def build_stress_trajectory_figure(
    records: list[Tier1aStressRunRecord],
    hedft: HedftTrajectory,
    smoothed: SmoothedSpeedReference,
    window: tuple[float, float],
    *,
    selected_n_final: int = PLOT_N_FINAL,
    title: str | None = None,
):
    import matplotlib.pyplot as plt

    fixed = [record for record in records if record.n_final is None]
    selected = [record for record in records if record.n_final == selected_n_final]
    if len(fixed) != 1 or len(selected) != 1:
        raise ValueError(
            f"expected one fixed and one stress n_final={selected_n_final} record"
        )

    fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    for record in [fixed[0], selected[0]]:
        t_md, _, _, v2_md = ensemble_mean_series(record.ion)
        ax.plot(t_md, v2_md, lw=1.4, ls="--", label=f"{record.label} MD mean |v2|")

    ax.plot(smoothed.time_ps, smoothed.speed_Aps, color="black", lw=1.6, ls=":", label="CEEMDAN+SG |v2|")
    ax.plot(hedft.time_ps, hedft.v2_magnitude_Aps, color="0.15", lw=1.2, alpha=0.65, label="HeDFT |v2|")
    ax.axvspan(window[0], window[1], color="tab:green", alpha=0.12, label="scored window")
    ax.set_ylabel(r"$|v|$ / $\mathrm{\AA}/\mathrm{ps}$")
    ax.set_xlabel("t / ps")
    ax.legend(frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle(title or f"Tier-1a onset-strip stress |v2| (n_final={selected_n_final})")
    return fig


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "(no Tier-1a stress rows)"
    widths = {
        col: max(len(col), *(len(_format_value(row[col])) for row in rows))
        for col in TIER1A_STRESS_TABLE_COLUMNS
    }
    header = "  ".join(col.ljust(widths[col]) for col in TIER1A_STRESS_TABLE_COLUMNS)
    sep = "  ".join("-" * widths[col] for col in TIER1A_STRESS_TABLE_COLUMNS)
    body = [
        "  ".join(_format_value(row[col]).ljust(widths[col]) for col in TIER1A_STRESS_TABLE_COLUMNS)
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIER1A_STRESS_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return p


def export_tier1a_stress_mean_series(
    directory: str | Path,
    records: list[Tier1aStressRunRecord],
) -> list[Path]:
    out_dir = Path(directory)
    paths: list[Path] = []
    for record in records:
        path = out_dir / f"{record.run_tag}_mean_trajectory.csv"
        provenance = (
            "Tier-1a onset-strip stress MD ensemble-mean trajectory. "
            f"run_tag={record.run_tag}, n_final={record.n_final}, "
            f"t_strip_ps={record.t_strip_ps}, N={record.ion.num_molecules}, "
            f"ion_steps={record.ion.time_ps.size}."
        )
        export_mean_series(path, record.ion, provenance=provenance)
        paths.append(path)
    return paths


def main() -> int:
    records, hedft, smoothed = collect_tier1a_stress_records()
    rows = [record.row for record in records if record.row]
    window = (WINDOW_START_PS, float(smoothed.time_ps[-1]))
    print("Tier-1a onset-strip stress table: diagnostic sensitivity runs only.")
    print(format_table(rows))
    if SAVE_CSV_PATH is not None:
        print(f"Wrote Tier-1a stress table -> {write_rows_csv(SAVE_CSV_PATH, rows)}")
    if EXPORT_MEAN_SERIES_DIR is not None:
        for path in export_tier1a_stress_mean_series(EXPORT_MEAN_SERIES_DIR, records):
            print(f"Wrote Tier-1a stress mean series -> {path}")
    if SHOW_FIGURE:
        import matplotlib.pyplot as plt

        build_stress_trajectory_figure(records, hedft, smoothed, window, selected_n_final=PLOT_N_FINAL)
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run script tests and compile**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier1a_scripts.py -q
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m py_compile scripts/post_processing/tier1a_stress_table.py
```

Expected: tests pass with the existing warning; py_compile exits 0.

- [ ] **Step 5: Commit Task 6**

```powershell
git add scripts/post_processing/tier1a_stress_table.py tests/test_tier1a_scripts.py
git commit -m "Add Tier-1a onset strip stress scorer"
```

---

### Task 7: Add Run-Level Stress Smoke

**Files:**
- Modify: `tests/test_ion_drag_smoke.py`

- [ ] **Step 1: Write full run smoke test**

In `tests/test_ion_drag_smoke.py`, add:

```python
def test_onset_strip_stress_runs_and_jumps_to_requested_shell():
    from dataclasses import replace
    from i2_helium_md.physics.shell_schedule import complex_mass_amu

    neutral = _synthetic_neutral(num_molecules=2, mass_amu=210.9546)
    cfg = replace(
        _anchored_cfg(),
        anchor_mode="onset_strip",
        anchor_n_final=0,
        t_star_ps=0.5,
        mass_initial_amu=complex_mass_amu(21),
    )
    ck = run_ion_propagation(cfg, neutral)

    n_shell = ck.n_shell[0, :]
    assert n_shell[0] == 21
    assert n_shell[-1] == 0
    assert np.array_equal(np.unique(n_shell), np.array([0, 21]))
    assert ck.E_mass_transfer_eV[:, -1].sum() > 0.0
    assert np.all(np.isfinite(ck.velocities_x))
```

- [ ] **Step 2: Run smoke tests**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_ion_drag_smoke.py -q
```

Expected: all smoke tests pass.

- [ ] **Step 3: Commit Task 7**

```powershell
git add tests/test_ion_drag_smoke.py
git commit -m "Add onset strip stress smoke coverage"
```

---

### Task 8: Update Docs

**Files:**
- Modify: `TIER1A_IMPLEMENTATION_PLAN.md`
- Modify: `drag_migration_log_tier1a.md`

- [ ] **Step 1: Update implementation plan**

In `TIER1A_IMPLEMENTATION_PLAN.md`, add a completed slice after Slice C:

```markdown
### Slice V -- Onset-violent stripping stress tests *(diagnostic; completed)*

> **IMPLEMENTED (2026-06-29).** Added a separate `onset_strip` stress schedule
> family with one violent event at `t=0.5 ps`, stripping from `n=21` to
> `n_final in {14,2,1,0}`. This is a sensitivity diagnostic only; it is not the
> TDDFT-anchored physical Tier-1a model.

The stress path keeps the corrected continuous-velocity shed:
`v_plus = v_minus`, with
`dE_mass_transfer = +0.5*(21-n_final)*m_He*|v|^2`.
Stress runs use separate tags, generator, scorer, and plots so they cannot be
confused with the physical `anchored_discrete` continuous baseline.
```

- [ ] **Step 2: Update migration log**

In `drag_migration_log_tier1a.md`, append:

```markdown
## Tier-1a Slice V -- delivery record (2026-06-29): onset-violent stripping stress tests

Implemented a separate diagnostic stress family, `onset_strip`, for violent
single-event shell stripping at `t=0.5 ps`. The endpoints are `n_final={14,2,1,0}`.
The mass update remains continuous-velocity and books the co-moving removed-He
kinetic energy positively in `E_mass_transfer_eV`.

These runs are sensitivity tests only. They do not supersede the physical
Tier-1a anchored continuous-velocity baseline and they do not reintroduce
cold-shed velocity kicks.
```

- [ ] **Step 3: Commit Task 8**

```powershell
git add TIER1A_IMPLEMENTATION_PLAN.md drag_migration_log_tier1a.md
git commit -m "Document Tier-1a onset stripping stress slice"
```

---

### Task 9: Final Verification

**Files:**
- No code edits unless verification exposes a bug.

- [ ] **Step 1: Run targeted test matrix**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_shell_schedule.py tests/test_mass_jump.py tests/test_ion_variable_mass.py tests/test_ion_drag_smoke.py tests/test_tier1a_scripts.py -q
```

Expected: all pass; the known `anchored_discrete` constant-coefficient pairing warning may appear.

- [ ] **Step 2: Run full suite**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

Expected: all tests pass; the known pairing warning may appear.

- [ ] **Step 3: Run static import checks**

Run:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m py_compile scripts/gen_tier1a_stress_runs.py scripts/post_processing/tier1a_stress_table.py
```

Expected: no output, exit code 0.

- [ ] **Step 4: Check stale naming**

Run:

```powershell
rg -n "tier1a_stress_onset_strip|onset_strip|build_onset_strip_schedule" i2_helium_md scripts tests TIER1A_IMPLEMENTATION_PLAN.md drag_migration_log_tier1a.md
```

Expected: references appear only in the new stress schedule, stress helpers/scripts, tests, and docs.

- [ ] **Step 5: Inspect final status**

Run:

```powershell
git status --short
```

Expected: only intentional project artifacts remain untracked or modified. If a
verification fix was required in this task, commit the exact files changed by that
fix with a specific message before handing off.
