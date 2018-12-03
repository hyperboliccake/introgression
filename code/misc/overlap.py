def overlap(start1, end1, start2, end2):
    if end1 < start2 or end2 < start1:
        return 0
    if start1 < start2:
        if start2 < end1:
            if end2 < end1:
                return end2 - start2 + 1
            else:
                return end1 - start2 + 1
        else:
            return 0
    else:
        if start1 < end2:
            if end1 < end2:
                return end1 - start1 + 1
            else:
                return end2 - start1 + 1
        else:
            return 0
        
    #if start1 < start2:
    #    return max(end1 - start2 + 1, 0) - max(end1 - end2, 0)
    #return max(end2 - start1 + 1, 0) - max(end2 - end1, 0)

def overlap_any(start1, end1, coords):
    for start2, end2 in coords:
        if overlap(start1, end1, start2, end2):
            return True
    return False

def contained(i, start, end):
    return i >= start and i <= end

def contained_any(i, coords):
    for start2, end2 in coords:
        if contained(i, start2, end2):
            return True
    return False

def contained_any_named(i, coords):
    for start2, end2 in coords.keys():
        if contained(i, start2, end2):
            return coords[(start2, end2)]
    return None

def overlap_region(start1, end1, start2, end2):
    o_start = max(start1, start2)
    o_end = min(end1, end2)
    if o_start > o_end:
        return -1, -1 # disjoint ranges
    return o_start, o_end
