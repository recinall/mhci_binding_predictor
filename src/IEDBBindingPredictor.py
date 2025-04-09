# This file is kept for backward compatibility
# The code has been moved to the mhci_binding_predictor package

from mhci_binding_predictor import IEDBBindingPredictor

# Import warning
import warnings
warnings.warn(
    "This module is deprecated. Please use 'from mhci_binding_predictor import IEDBBindingPredictor' instead.",
    DeprecationWarning,
    stacklevel=2
)
