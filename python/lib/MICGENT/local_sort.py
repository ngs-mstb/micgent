"""Sort query by reference when reodering is confined to a window"""

def sort_by_ref_in_window(inp_ref,inp_que,window_size=1000):
    """Reorder inp_que iterator by inp_ref within window_size.
    
    Assumptions:
    - all elements of inp_que are in inp_ref
    - order of inp_que might diverge from order of inp_ref
      only within window_size positions of inp_que

    Inputs:
    - inp_ref - iterator of IDs
    - inp_que - iterator of (ID, body) tuples for each element

    Yield:
    - (ID, body) tuples from inp_que

    Complexity:
    - Time is linear in both inp_que and inp_ref
    - Memory is capped by window_size inp_que elements

    Exceptions:
    - Raise ValueError with a dict containing any remaining
      elements of inp_queue (capped by window_size) after
      inp_ref is exhausted

    Applications:
    - Provide determinism to multi-threaded sequence filters /
      transformers such as bbduk. Results from different threads
      can be mixed with slightly different order, but the
      differences should be localized within a certain window
      around each position in the original input sequence ordering.
    """
    d_que = dict()
    ## prefetch windows_size of query
    for i_que in range(window_size):
        try:
            el_id, el_body = next(inp_que)
            d_que[el_id] = el_body
        except StopIteration:
            break
    ## stream through entire reference, output every
    ## matching query element from the cache and immediately 
    ## replace the matched cache element with a newly fetched query element
    for el_id in inp_ref:
        el_body = d_que.get(el_id,None)
        if el_body is not None:
            print(el_id)
            yield (el_id,el_body)
            del d_que[el_id]
            try:
                el_id, el_body = next(inp_que)
                d_que[el_id] = el_body
            except StopIteration:
                pass
    ## if anything is still incoming from inp_que,
    ## fetch it up to window_size for error reporting
    for i_que in range(window_size-len(d_que)):
        try:
            el_id, el_body = next(inp_que)
            d_que[el_id] = el_body
        except StopIteration:
            break
    if len(d_que) > 0:
        raise ValueError(d_que)

def test_sort_by_ref_in_window():
    inp_ref = [1,2,3,4,5,6,7,8,9]
    inp_que = [2,3,6,4,7,8,9]
    out_exp = [2,3,4,6,7,8,9]
    window_size = 2
    out_obs = []
    for el_id, el_body in sort_by_ref_in_window(inp_ref=inp_ref,
        inp_que=zip(inp_que,inp_que),
        window_size=window_size):
        out_obs.append(el_id)
        #print([el_id,el_body])
    assert(out_obs == out_exp)

