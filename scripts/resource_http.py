"""Small REST client for provenance-preserving resource downloads."""
from pathlib import Path
import json

import requests


def execute(payload):
    url = payload['base_url'].rstrip('/') + '/' + payload['path'].lstrip('/')
    try:
        with requests.get(url, params=payload.get('params'),
                          timeout=payload.get('timeout_sec', 180)) as response:
            response.raise_for_status()
            if payload.get('response_format') == 'text':
                raw = response.text
            else:
                raw = json.dumps(response.json(), indent=2)
            path = Path(payload['raw_output_path'])
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(raw, encoding='utf-8')
            return {'ok': True, 'raw_output_path': str(path),
                    'status_code': response.status_code}
    except (requests.RequestException,ValueError) as exc:
        return {'ok': False, 'error': str(exc)}
