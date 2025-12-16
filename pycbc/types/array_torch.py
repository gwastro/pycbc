# array_torch.py

import torch
import numpy as np
import pycbc.scheme as _scheme
from pycbc.types.array import Array

# Mapping from NumPy dtypes to Torch dtypes
NUMPY_TO_TORCH_DTYPE = {
    np.dtype(np.float32): torch.float32,
    np.dtype(np.float64): torch.float64,
    np.dtype(np.complex64): torch.complex64,
    np.dtype(np.complex128): torch.complex128,
    np.dtype(np.int32): torch.int32,
    np.dtype(np.uint32): torch.uint32,
    np.dtype(int): torch.int64,
}

# Reverse mapping from Torch dtypes to NumPy dtypes
TORCH_TO_NUMPY_DTYPE = {v: k for k, v in NUMPY_TO_TORCH_DTYPE.items()}

def _scheme_get_numpy_dtype(dtype):
    """
    Get the NumPy dtype corresponding to the given Torch dtype.

    Parameters
    ----------
    dtype : torch.dtype
        The Torch dtype to convert.

    Returns
    -------
    numpy.dtype
        The corresponding NumPy dtype.
    """
    if isinstance(dtype, torch.dtype):
        numpy_dtype = TORCH_TO_NUMPY_DTYPE.get(dtype)
        if numpy_dtype is None:
            raise TypeError(f"Torch data type {dtype} does not have a corresponding NumPy dtype")
        return numpy_dtype
    else:
        # If it's already a NumPy dtype, return it
        return np.dtype(dtype)

def _scheme_array_from_initial(initial_array, dtype=None):
    device = _scheme.mgr.state.device

    if isinstance(initial_array, torch.Tensor):
        tensor = initial_array.to(device)
        if dtype is not None:
            torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(dtype))
            if tensor.dtype != torch_dtype:
                tensor = tensor.to(dtype=torch_dtype)
    else:
        # Convert initial_array to torch.Tensor
        if not isinstance(initial_array, np.ndarray):
            initial_array = np.array(initial_array, dtype=dtype)
        else:
            if dtype is not None and initial_array.dtype != dtype:
                initial_array = initial_array.astype(dtype)
        torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(initial_array.dtype))
        tensor = torch.tensor(initial_array, dtype=torch_dtype, device=device)

    return tensor

def _scheme_matches_base_array(array):
    return isinstance(array, torch.Tensor)

def _copy_base_array(array):
    return array.clone()

def _to_device(array):
    """
    Move the array to the appropriate device and convert it to a Torch tensor.

    Parameters
    ----------
    array : array-like
        Input data, can be a NumPy array, list, or Torch tensor.

    Returns
    -------
    torch.Tensor
        Torch tensor on the specified device.
    """
    device = _scheme.mgr.state.device

    if isinstance(array, torch.Tensor):
        tensor = array.to(device)
    else:
        # Convert array to Torch tensor
        tensor = torch.tensor(np.asanyarray(array))
        tensor = tensor.to(device)
    return tensor

def zeros(length, dtype=np.float64):
    torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(dtype))
    device = _scheme.mgr.state.device
    return Array(torch.zeros(length, dtype=torch_dtype, device=device), copy=False)

def empty(length, dtype=np.float64):
    torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(dtype))
    device = _scheme.mgr.state.device
    return Array(torch.empty(length, dtype=torch_dtype, device=device), copy=False)

def ptr(array):
    return array._data.data_ptr()

def numpy(array):
    """
    Convert the Array to a NumPy array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    numpy.ndarray
        NumPy array.
    """
    return array._data.cpu().numpy()

def inner(a, b):
    """
    Compute the inner product of two arrays.

    Parameters
    ----------
    a, b : Array
        Input PyCBC Array instances.

    Returns
    -------
    scalar
        The inner product result.
    """
    data_a = a._data
    data_b = b
    if data_a.is_complex():
        result = torch.vdot(data_a.view(-1), data_b.view(-1))
    else:
        result = torch.dot(data_a.view(-1), data_b.view(-1))
    return result.item()

vdot = inner  # Alias for inner product

def squared_norm(array):
    """
    Compute the squared norm (magnitude squared) of the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    Array
        An Array instance containing the squared norm.
    """
    data = array._data
    result_data = torch.abs(data) ** 2
    return Array(result_data, copy=False)

def sum(array):
    """
    Compute the sum of all elements in the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    scalar
        Sum of elements.
    """
    return torch.sum(array._data).item()

def max(array):
    """
    Find the maximum value in the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    scalar
        Maximum value.
    """
    return torch.max(array._data).item()

def max_loc(array):
    """
    Find the maximum value and its index in the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    tuple
        Maximum value and its index.
    """
    data = array._data
    max_val, max_idx = torch.max(data, dim=0)
    return max_val.item(), max_idx.item()

def min(array):
    """
    Find the minimum value in the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    scalar
        Minimum value.
    """
    return torch.min(array._data).item()

def cumsum(array):
    """
    Compute the cumulative sum of the array along dimension 0.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    Array
        An Array instance containing the cumulative sum.
    """
    data = array._data
    result_data = torch.cumsum(data, dim=0)
    return Array(result_data, copy=False)

def take(array, indices):
    """
    Extract elements from the array at the specified indices.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    indices : array-like
        Indices to take.

    Returns
    -------
    Array
        An Array instance containing the extracted elements.
    """
    data = array._data
    indices_tensor = torch.tensor(indices, dtype=torch.long, device=data.device)
    result_data = torch.take(data, indices_tensor)
    return Array(result_data, copy=False)

def dot(a, b):
    """
    Compute the dot product of two arrays.

    Parameters
    ----------
    a, b : Array
        Input PyCBC Array instances.

    Returns
    -------
    scalar
        Dot product result.
    """
    data_a = a._data
    data_b = b
    result = torch.dot(data_a.view(-1), data_b.view(-1))
    return result.item()

def abs_max_loc(array):
    """
    Find the maximum absolute value and its index in the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    tuple
        Maximum absolute value and its index.
    """
    data = array._data
    abs_data = torch.abs(data)
    max_val, max_idx = torch.max(abs_data, dim=0)
    return max_val.item(), max_idx.item()

def clear(array):
    """
    Zero out all elements of the array in-place.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    None
    """
    array._data.zero_()

def _getvalue(array, index):
    """
    Get the value at the specified index from the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    index : int
        Index to retrieve.

    Returns
    -------
    scalar
        Value at the specified index.
    """
    return array._data[index].item()

def _copy(dest, src):
    """
    Copy data from src array to dest array.

    Parameters
    ----------
    dest, src : Array
        Destination and source PyCBC Array instances.

    Returns
    -------
    None
    """
    dest._data.copy_(src._data)

def multiply_and_add(a, b, mult_fac):
    """
    Multiply b by mult_fac and add to a.

    Parameters
    ----------
    a, b : Array
        Input PyCBC Array instances.
    mult_fac : scalar
        Multiplication factor.

    Returns
    -------
    None
    """
    a._data.add_(b._data, alpha=mult_fac)

def weighted_inner(a, b, weight):
    """
    Compute the weighted inner product of two arrays.

    Parameters
    ----------
    a, b : Array
        Input PyCBC Array instances.
    weight : Array or None
        Weights array or None.

    Returns
    -------
    scalar
        Weighted inner product result.
    """
    data_a = a._data
    data_b = b
    if weight is None:
        return inner(a, b)
    else:
        data_w = weight._data
        if data_a.is_complex():
            result = torch.dot((data_a.conj() * data_b / data_w).view(-1), torch.ones_like(data_a.view(-1)))
        else:
            result = torch.dot((data_a * data_b / data_w).view(-1), torch.ones_like(data_a.view(-1)))
        return result.item()

def real(array):
    """
    Return the real part of the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    Array
        An Array instance containing the real part.
    """
    data = torch.real(array._data)
    return Array(data, copy=False)

def imag(array):
    """
    Return the imaginary part of the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    Array
        An Array instance containing the imaginary part.
    """
    data = torch.imag(array._data)
    return Array(data, copy=False)

def conj(array):
    """
    Return the complex conjugate of the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.

    Returns
    -------
    Array
        An Array instance containing the complex conjugate.
    """
    data = torch.conj(array._data)
    return Array(data, copy=False)

def fill(array, value):
    """
    Fill the array with the specified value.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    value : scalar
        Value to fill the array with.

    Returns
    -------
    None
    """
    array._data.fill_(value)

def roll(array, shift):
    """
    Roll the elements of the array by the specified shift.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    shift : int
        Number of places by which elements are shifted.

    Returns
    -------
    Array
        An Array instance with rolled elements.
    """
    data = torch.roll(array._data, shifts=shift, dims=0)
    return Array(data, copy=False)

def astype(array, dtype):
    """
    Change the data type of the array.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    dtype : numpy.dtype
        Desired NumPy data type.

    Returns
    -------
    Array
        An Array instance with the new data type.
    """
    torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(dtype))
    data = array._data.to(dtype=torch_dtype)
    return Array(data, copy=False)

def resize(array, new_size):
    """
    Resize the array to the specified new size.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    new_size : int
        New size of the array.

    Returns
    -------
    Array
        An Array instance with the new size.
    """
    data = array._data
    current_size = data.size(0)
    if new_size == current_size:
        return array
    elif new_size < current_size:
        new_data = data[:new_size].clone()
    else:
        # Create a new tensor with the new size and copy existing data
        torch_dtype = data.dtype
        device = data.device
        new_data = torch.empty(new_size, dtype=torch_dtype, device=device)
        new_data[:current_size] = data
    return Array(new_data, copy=False)

def view(array, dtype):
    """
    Return a view of the array with the specified data type.

    Parameters
    ----------
    array : Array
        Input PyCBC Array instance.
    dtype : numpy.dtype
        Desired NumPy data type.

    Returns
    -------
    Array
        An Array instance with the new data type.
    """
    torch_dtype = NUMPY_TO_TORCH_DTYPE.get(np.dtype(dtype))
    data = array._data.to(dtype=torch_dtype)
    return Array(data, copy=False)

# Any other functions you need can be added following the same pattern.

