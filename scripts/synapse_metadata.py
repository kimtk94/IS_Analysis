"""Bounded-concurrency Synapse metadata helpers for user-run batch workflows."""
from __future__ import annotations

import asyncio
import importlib.util
import os
from typing import Any, Callable


def fetch_file_metadata_many(
    synapse_ids: list[str], token: str = "", max_concurrency: int = 8,
    progress: Callable[[int, int], None] | None = None,
) -> dict[str, dict[str, str]]:
    """Fetch size and checksums only for selected files, with bounded concurrency."""
    if importlib.util.find_spec("synapseclient") is None:
        raise RuntimeError("synapseclient is required for Synapse metadata lookup.")
    if max_concurrency < 1:
        raise ValueError("max_concurrency must be at least 1")
    import synapseclient
    from synapseclient.api import get_entity
    from synapseclient.api.file_services import get_file_handle_for_download_async

    unique_ids = list(dict.fromkeys(item for item in synapse_ids if item))
    syn = synapseclient.Synapse()
    syn.login(authToken=token or None, silent=not bool(token))

    async def collect() -> dict[str, dict[str, str]]:
        semaphore = asyncio.Semaphore(max_concurrency)
        completed = 0

        async def fetch(synapse_id: str) -> tuple[str, dict[str, str]]:
            nonlocal completed
            async with semaphore:
                entity = await get_entity(synapse_id, synapse_client=syn)
                handle_id = entity.get("dataFileHandleId")
                if not handle_id:
                    raise RuntimeError(f"Synapse entity has no dataFileHandleId: {synapse_id}")
                info = await get_file_handle_for_download_async(
                    str(handle_id), synapse_id, synapse_client=syn,
                )
                handle = info.get("fileHandle", {})
            completed += 1
            if progress and (completed % 100 == 0 or completed == len(unique_ids)):
                progress(completed, len(unique_ids))
            return synapse_id, {
                "expected_size_bytes": str(handle.get("contentSize", "")),
                "md5": handle.get("contentMd5", ""),
                "sha256": handle.get("contentSha256", ""),
            }

        pairs = await asyncio.gather(*(fetch(synapse_id) for synapse_id in unique_ids))
        return dict(pairs)

    return asyncio.run(collect())


def synapse_token() -> str:
    """Use the token convention shared by the manifest builder and runner."""
    return os.environ.get("SYNAPSE_AUTH_TOKEN", "")
