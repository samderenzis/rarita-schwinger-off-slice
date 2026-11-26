"""Run ONE unmodified parent check script with rs_model standing in for the parent's bhmap_defs.

Usage (called by step1_slice_regression.py, one fresh process per script):
    python3 -I parent_shim_runner.py <path/to/parent/check_N_*.py>

The stand-in module is rs_model, except that symbolic_model() and witness_model() pass the charge
EXPLICITLY, q = lam*sqrt(2/3)*m, computed here from the model's own m and lam. This exercises the new
"q given as a parameter" code path, not the q=None default. The parent script's own folder is NOT put
on sys.path, so the parent bhmap_defs.py cannot be imported by accident.
"""
import os
import runpy
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp  # noqa: E402
import rs_model  # noqa: E402


def explicit_slice_q(m, lam):
    # written out independently of rs_model.pauli_q on purpose
    return lam * sp.sqrt(sp.Rational(2, 3)) * m


def symbolic_model():
    m, lam, b1, b2, pip = sp.symbols('m lam b1 b2 pip', positive=True)
    M = rs_model.Model(m, lam, b1, b2, pip, q=explicit_slice_q(m, lam))
    assert not M.on_slice  # q was passed explicitly
    return M


def witness_model():
    m, lam = sp.sqrt(sp.Rational(823, 142944)), sp.sqrt(sp.Rational(1, 2))
    M = rs_model.Model(m, lam, 1, sp.Rational(1, 3), sp.Rational(1, 2), q=explicit_slice_q(m, lam))
    assert not M.on_slice
    return M


shim = types.ModuleType('bhmap_defs')
shim.__dict__.update({k: v for k, v in vars(rs_model).items() if not k.startswith('__')})
shim.symbolic_model = symbolic_model
shim.witness_model = witness_model
sys.modules['bhmap_defs'] = shim

target = os.path.abspath(sys.argv[1])
sys.argv = [target]
runpy.run_path(target, run_name='__main__')
