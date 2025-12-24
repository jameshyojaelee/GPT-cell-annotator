from __future__ import annotations

from backend.llm import prompts


def test_single_cluster_prompt_includes_markers():
    cluster = {
        "cluster_id": "7",
        "markers": ["MS4A1"],
    }
    prompt = prompts.build_single_cluster_prompt(cluster, {"species": "Homo sapiens"})
    assert "Cluster ID: 7" in prompt
    assert "MS4A1" in prompt


def test_batch_prompt_formats_clusters():
    clusters = [
        {
            "cluster_id": "0",
            "markers": ["CD3D"],
        }
    ]
    prompt = prompts.build_batch_prompt(clusters, {"species": "Homo sapiens"})
    assert "Cluster 0" in prompt
    assert "CD3D" in prompt
