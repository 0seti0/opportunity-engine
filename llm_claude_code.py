"""LLM judgment via Claude Code headless (`claude -p`) — runs on the Max SUBSCRIPTION, NO API key.

Shells out to the `claude` CLI (which is authenticated by your subscription), tells the model to return
schema-shaped JSON, and parses it with retries. Swap this one file for an API-key client later (same
`judge()` signature) when you productize for customers.
"""
import json
import subprocess


def _extract_json(text):
    t = (text or "").strip()
    if t.startswith("```"):                       # strip a ```json fence if the model added one
        t = t.strip("`")
        t = t[t.find("\n") + 1:] if "\n" in t else t
    candidates = [t]
    if "{" in t and "}" in t:
        candidates.append(t[t.find("{"): t.rfind("}") + 1])
    for c in candidates:
        try:
            return json.loads(c)
        except Exception:
            pass
    return None


def judge(evidence: dict, schema: dict, model: str | None = None, retries: int = 2) -> dict:
    """Reason over a pre-assembled evidence bundle -> schema-valid verdict. No web search, no facts to
    look up (those are in `evidence`); the model only WEIGHS the evidence."""
    keys = ", ".join(schema.get("properties", {}).keys())
    prompt = (
        "You are a covalent drug-discovery analyst. Judge this target using ONLY the EVIDENCE below. "
        "Do NOT restate residue or selectivity facts — those are verified deterministically and will "
        "overwrite anything you say about them. Your job is to WEIGH the evidence into a verdict.\n\n"
        f"Return ONLY one JSON object with exactly these keys: {keys}. No markdown, no prose, no fences.\n\n"
        f"SCHEMA: {json.dumps(schema)}\n\nEVIDENCE: {json.dumps(evidence)}"
    )
    cmd = ["claude", "-p", prompt, "--output-format", "json"]
    if model:
        cmd += ["--model", model]
    last = None
    for _ in range(retries + 1):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            env = json.loads(r.stdout)            # `claude -p --output-format json` envelope
            out = _extract_json(env.get("result", ""))
            if out and all(k in out for k in schema.get("required", [])):
                out["_cost_usd"] = env.get("total_cost_usd", 0)   # 0 on subscription
                return out
            last = env.get("result", r.stdout)
        except Exception as e:
            last = str(e)
    raise RuntimeError(f"claude -p judgment failed: {str(last)[:300]}")


if __name__ == "__main__":                        # smoke test the no-key plumbing
    s = {"type": "object", "properties": {"ok": {"type": "boolean"}, "msg": {"type": "string"}},
         "required": ["ok", "msg"]}
    print(judge({"ping": "hello"}, s))
