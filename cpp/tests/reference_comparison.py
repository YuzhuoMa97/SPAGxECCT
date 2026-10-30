"""Separate CGF-dependent outputs from unchanged model and score quantities."""
import math


def depends_on_spa_cgf(name):
    return name.startswith(("LOG10P_Gx", "LOG10P_SPA_Gx", "Z_Gx",
                            "SPA_STATUS_Gx", "LOG10P_ANCESTRY_CCT_Gx"))


def comparison(current, reference, tolerance):
    if current == reference:
        return True, 0.0
    try:
        x, y = float(current), float(reference)
    except ValueError:
        return False, None
    if math.isnan(x) and math.isnan(y):
        return True, 0.0
    if not math.isfinite(x) or not math.isfinite(y):
        return False, None
    error = abs(x-y)
    return error <= tolerance, error


def valid_changed_p(name, current, reference):
    """A legal approximation change must not silently lose a formerly valid P."""
    if not name.startswith("LOG10P_"):
        return True
    try:
        old = float(reference)
    except ValueError:
        return True
    if math.isnan(old):
        return True
    try:
        value = float(current)
    except ValueError:
        return False
    return not math.isnan(value) and value >= 0
