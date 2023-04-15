"""
Author: Yun Zhang, yzhxdy@gmail.com
Update: 2018

substring search using Rabin-Karp algorithm
text: document to be searched
pattern: keyword, length >= 1
      if len(pattern) is 0, hash_value = 0, precomputed hash values for all substrings (hash_substrings) has 1 more position.
prime: a very large number to minimize collision. Take everything mod prime as soon as possible to keep the number < prime.
      Be aware of integer overflow.
"""
import random


def polyhash(string, prime, x):
    '''
    :return: hash value of input string
    '''
    hash_value = 0
    for char in reversed(string):
        hash_value = (hash_value * x + ord(char)) % prime
    return hash_value


def precompute_hashes(text, len_p, prime, x):
    '''
    uses modular hashing to compute hash values in linear time
    :param len_p: >= 1, otherwise, resulting array has one extra position
    :param prime: a very large number to decrease probability of collision
    :param x: 1...prime-1
    :return: an array of hash values for each substring of len_p.
    '''
    #x = random.randint(1, prime)
    len_t = len(text)
    # constraints: len_p >= 1, otherwise, this is 1+ off
    hash_substrings = [None] * (len_t - len_p + 1)
    last_substring = text[(len_t - len_p):]
    hash_substrings[len_t - len_p] = polyhash(last_substring, prime, x)

    y = 1
    for i in range(1, len_p + 1): # 1 to len_p
        y = (y * x) % prime

    for j in range((len_t - len_p - 1), -1, -1): # len_t - len_p - 1 to 0
        hash_substrings[j] = (x * hash_substrings[j + 1] + ord(text[j]) - y * ord(text[j + len_p])) % prime
    return hash_substrings


def rabin_karp_substring_search(pattern, text):
    prime = 1000000007 # prime: a very large number to decrease probability of collision
    x = random.randint(1, prime) # 1...prime-1
    #x = 263
    len_p = len(pattern)
    len_t = len(text)
    indices = []

    hash_pattern = polyhash(pattern, prime, x)
    # print('hash_pattern =', hash_pattern)
    # precompute hashes for all substrings of len_p, O(n)
    hash_substrings = precompute_hashes(text, len_p, prime, x)
    # print('hash_substrings =', hash_substrings)
    # scan all substrings of len_p
    for i in range(len_t - len_p + 1):
        if hash_pattern == hash_substrings[i]:
            if text[i:i + len_p] == pattern:
                indices.append(i)
    return indices


if __name__ == '__main__':
    pattern = 'abc'
    text = 'abcdbdabc'
    print(rabin_karp_substring_search(pattern, text))
