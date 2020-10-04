def rv_ran4_generator(seed):
    i4_huge = 2147483647
    seed_num = seed
    seed_num = round(seed,0)
    seed_num = seed_num % 2147483647

    if seed_num < 0:
        seed_num = seed_num + i4_huge

    k = round(seed_num/127773,0)
    seed_num = 16807 * (seed_num - k * 127773) - k * 2836

    if seed_num < 0:
        seed_num = seed_num + i4_huge

    value = seed_num * 4.656612875E-10
    return value


def rv_ran4_seed(seed):
    i4_huge = 2147483647
    seed_num = seed
    seed_num = round(seed,0)
    seed_num = seed_num % 2147483647

    if seed_num < 0:
        seed_num = seed_num + i4_huge

    k = round(seed_num/127773,0)
    seed_num = 16807 * (seed_num - k * 127773) - k * 2836

    if seed_num < 0:
        seed_num = seed_num + i4_huge

    value = seed_num
    return value
