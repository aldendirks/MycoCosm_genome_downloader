import os
from pathlib import Path
import subprocess as sp

import pytest # type: ignore

from mycocosm_downloader.auth import JGI_login


# ---------------------------------------------------
# helpers
# ---------------------------------------------------

class DummyProc:
    def check_returncode(self):
        return None


def fake_run_success(cmd, stdout=None, stderr=None):
    # simulate curl creating cookie file
    cookie_index = cmd.index("-c") + 1
    cookie_path = Path(cmd[cookie_index])
    cookie_path.write_text("cookie")
    return DummyProc()


def fake_run_fail(*args, **kwargs):
    raise sp.CalledProcessError(returncode=1, cmd="curl")


# ---------------------------------------------------
# tests
# ---------------------------------------------------

def test_login_with_env_vars(monkeypatch, tmp_path):
    cookie = tmp_path / "cookie.txt"

    monkeypatch.setenv("JGI_USER", "user1")
    monkeypatch.setenv("JGI_PASSWORD", "pass1")
    monkeypatch.setattr(sp, "run", fake_run_success)

    ok = JGI_login(cookie)

    assert ok is True
    assert cookie.exists()


def test_login_with_credentials_file(monkeypatch, tmp_path):
    cookie = tmp_path / "cookie.txt"
    creds = tmp_path / "creds.txt"
    creds.write_text("user2\npass2\n")

    monkeypatch.delenv("JGI_USER", raising=False)
    monkeypatch.delenv("JGI_PASSWORD", raising=False)
    monkeypatch.setattr(sp, "run", fake_run_success)

    ok = JGI_login(cookie, credentials_file=creds)

    assert ok is True
    assert cookie.exists()


def test_login_interactive(monkeypatch, tmp_path):
    cookie = tmp_path / "cookie.txt"

    monkeypatch.delenv("JGI_USER", raising=False)
    monkeypatch.delenv("JGI_PASSWORD", raising=False)

    inputs = iter(["user3", "n"])  # username, save_choice

    monkeypatch.setattr("builtins.input", lambda _: next(inputs))
    monkeypatch.setattr("mycocosm_downloader.auth.getpass", lambda _: "pass3")
    monkeypatch.setattr(sp, "run", fake_run_success)

    ok = JGI_login(cookie)

    assert ok is True
    assert cookie.exists()


def test_curl_failure(monkeypatch, tmp_path):
    cookie = tmp_path / "cookie.txt"

    monkeypatch.setenv("JGI_USER", "user")
    monkeypatch.setenv("JGI_PASSWORD", "pass")
    monkeypatch.setattr(sp, "run", fake_run_fail)

    ok = JGI_login(cookie)

    assert ok is False
