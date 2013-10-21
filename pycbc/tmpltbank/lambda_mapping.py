def generate_mapping(order):
    """
    This function will take an order string and return a mapping between
    components in the metric and the various Lambda components. This must be
    used (and consistently used) when generating the metric *and* when
    transforming to/from the xi_i coordinates to the lambda_i coordinates.

    NOTE: This is not a great way of doing this. It would be nice to clean
    this up. Hence pulling this function out.

    Parameters
    ----------
    order : string
        A string containing a PN order. Valid values are
        'zeroPN'
        'onePN'
        'onePointFivePN'
        'twoPN'
        'twoPointFivePN'
        'threePN'
        'threePointFivePN'
        'taylorF4_45PN'
    Returns
    --------
    mapping: dictionary
        A mapping between the active Lambda terms and index in the metric
    """
    mapping = {}
    if order == 'taylorF4_45PN':
        mapping['Lambda0'] = 0
        mapping['Lambda2'] = 1
        mapping['Lambda3'] = 2
        mapping['Lambda4'] = 3
        mapping['LogLambda5'] = 4
        mapping['Lambda6'] = 5
        mapping['LogLambda6'] = 6
        mapping['Lambda7'] = 7
        mapping['LogLambda8'] = 8
        mapping['LogLogLambda8'] = 9
        mapping['Lambda9'] = 10
        mapping['LogLambda9'] = 11
        return mapping
    else:
        mapping['Lambda0'] = 0
        if order == 'zeroPN':
            return mapping
        mapping['Lambda2'] = 1
        if order == 'onePN':
            return mapping
        mapping['Lambda3'] = 2
        if order == 'onePointFivePN':
            return mapping
        mapping['Lambda4'] = 3
        if order == 'twoPN':
            return mapping
        mapping['LogLambda5'] = 4
        if order == 'twoPointFivePN':
            return mapping
        mapping['Lambda6'] = 5
        mapping['LogLambda6'] = 6
        if order == 'threePN':
            return mapping
        mapping['Lambda7'] = 7
        if order == 'threePointFivePN':
            return mapping
        raise ValueError("Order %s is not understood." %(order))
