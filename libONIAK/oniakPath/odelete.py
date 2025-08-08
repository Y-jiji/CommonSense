import os, re, sys

""" Delete all files under the pattern. Returns number of files deleted. """
def delete_pattern(directory, pattern):
    cnt = 0
    for root, dirs, files in os.walk(directory):
        for file in files:
            full_path = os.path.join(root, file)
            if re.search(pattern, full_path):
                os.remove(full_path)
                cnt += 1
    return cnt

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python odelete.py [directory] [pattern]")
        exit(1)
    directory = sys.argv[1]
    pattern = sys.argv[2]
    cnt = delete_pattern(directory, pattern)
    print("Deleted {} files under pattern {}".format(cnt, pattern))