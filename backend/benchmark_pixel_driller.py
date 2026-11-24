#!/usr/bin/env python3
"""Benchmark pixel driller API

Usage:
    $ python benchmark_pixel_driller.py -c 4 -n 100
    Requests: 100 | Concurrency: 4
    Mean latency: 1543.54 ms
    Total wall time (approx): 154.35 s
    P95 latency: 2276.02 ms
    P99 latency: 2416.50 ms
    Status counts:
    200: 100
"""
import argparse
import asyncio
import random
import statistics
import time
from collections import Counter

import aiohttp


async def hit(
    session: aiohttp.ClientSession, base: str, sem: asyncio.Semaphore
) -> float:
    async with sem:
        lat = random.uniform(-90, 90)
        lon = random.uniform(-180, 180)
        url = f"{base.rstrip('/')}/{lon:.6f}/{lat:.6f}"
        start = time.perf_counter()
        async with session.get(url) as resp:
            await resp.read()
            elapsed = time.perf_counter() - start
            print(resp.status, url)
            return resp.status, elapsed


async def run(total: int, concurrency: int, base_url: str, timeout: float) -> None:
    sem = asyncio.Semaphore(concurrency)
    timeout_cfg = aiohttp.ClientTimeout(
        total=None, sock_connect=timeout, sock_read=timeout
    )
    async with aiohttp.ClientSession(timeout=timeout_cfg) as session:
        tasks = [asyncio.create_task(hit(session, base_url, sem)) for _ in range(total)]
        timings, statuses = [], Counter()
        for coro in asyncio.as_completed(tasks):
            status, elapsed = await coro
            statuses[status] += 1
            timings.append(elapsed)

    total_time = sum(timings)
    mean_latency = statistics.mean(timings)

    print(f"Requests: {total} | Concurrency: {concurrency}")
    print(f"Mean latency: {mean_latency*1000:.2f} ms")
    print(f"Total wall time (approx): {total_time:.2f} s")

    if total >= 20:
        p95 = statistics.quantiles(timings, n=100)[94]
        print(f"P95 latency: {p95*1000:.2f} ms")
    if total >= 100:
        p99 = statistics.quantiles(timings, n=100)[98]
        print(f"P99 latency: {p99*1000:.2f} ms")

    print("Status counts:")
    for code, count in sorted(statuses.items()):
        print(f"  {code}: {count}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark for pixel-driller endpoint."
    )
    parser.add_argument(
        "--base-url",
        default="http://localhost:8888/pixel-driller/point",
        help="Base URL up to /point (/lon/lat added with random values).",
    )
    parser.add_argument(
        "-n",
        "--requests",
        type=int,
        default=200,
        help="Total number of requests.",
    )
    parser.add_argument(
        "-c",
        "--concurrency",
        type=int,
        default=16,
        help="Number of concurrent requests.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=15.0,
        help="Per-request timeout in seconds.",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Seed for random coordinates."
    )
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    asyncio.run(run(args.requests, args.concurrency, args.base_url, args.timeout))


if __name__ == "__main__":
    main()
