"""Performs a grid search over two hyperparameters."""

next_step = {10:5, 5:3, 4:2, 3:2, 2:1, 1:1}

def range_min(ran, val):
    return max(ran[0], val)
def range_max(ran, val):
    return min(ran[1], val)

def grid_search(xrange, yrange, func, feedback, step=10, table={}):
    flag = True
    while flag:
        # grid run
        if step == 1:
            flag = False
        for x in range(*xrange, step):
            for y in range(*yrange, step):
                if (x, y) not in table.keys():
                    result = feedback(func(x, y))
                    if result is not None:
                        flag = True
                        print(x, y, result)
                        table[x, y] = result

        # find minimum 
        x, y = min(table, key=table.get)
        # if x == xrange[0] or x == xrange[1]-1 or y == yrange[0] or y == yrange[1]-1:
        #     print("Grid search aborted. Minimum found at boundary.")
        #     return table

        step = next_step[step]
        xrange = (range_min(xrange, x - 2 * step), range_max(xrange, x + 2 * step+1))
        yrange = (range_min(yrange, y - 2 * step), range_max(yrange, y + 2 * step+1))
    return x, y, table