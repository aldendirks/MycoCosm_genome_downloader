import os
import subprocess as sp
from getpass import getpass
from pathlib import Path
from platformdirs import user_config_dir
from subprocess import DEVNULL


def JGI_login(cookie_path: Path, credentials_file: Path | None = None) -> bool:
    """Login to JGI and create cookie file."""
    
    user = None
    password = None

    env_user = os.getenv("JGI_USER")
    env_pass = os.getenv("JGI_PASSWORD")

    if env_user and env_pass:
        user = env_user.strip()
        password = env_pass.strip()
        print("Using JGI credentials from environment variables.")

    if (not user or not password) and credentials_file and credentials_file.exists():
        try:
            lines = credentials_file.read_text().splitlines()
            if len(lines) < 2:
                raise ValueError
            user, password = lines[0].strip(), lines[1].strip()
            print(f"Using credentials from {credentials_file}.")
        except Exception:
            print("Invalid credentials file — will prompt manually.")

    if not user or not password:
        print("Please enter your JGI credentials:")
        user = input("Username: ").strip()
        password = getpass("Password: ").strip()

        if not user or not password:
            raise RuntimeError("Missing JGI credentials.")
        
        save_choice = input("\nSave credentials to config file for future use? (yes/no): ").strip().lower()

        if save_choice in ("yes", "y"):
            config_dir = Path(user_config_dir("mycocosm_downloader"))
            config_dir.mkdir(parents=True, exist_ok=True)

            creds_file = config_dir / "credentials.txt"
            creds_file.write_text(f"{user}\n{password}\n")

            try:
                creds_file.chmod(0o600)
            except Exception:
                pass

            print(f"Credentials saved to {creds_file}.")

            print("\nTip: You can instead use environment variables:")
            print("  export JGI_USER=your_username")
            print("  export JGI_PASSWORD=your_password")

    cookie_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "curl",
        "--silent",
        "https://signon.jgi.doe.gov/signon/create",
        "--data-urlencode", f"login={user}",
        "--data-urlencode", f"password={password}",
        "-c", str(cookie_path),
    ]

    try:
        proc = sp.run(cmd, stdout=DEVNULL, stderr=DEVNULL)
        proc.check_returncode()
    except FileNotFoundError:
        print("Error: curl is not installed or not in PATH.")
        return False
    except sp.CalledProcessError:
        print("JGI login request failed.")
        return False

    if not cookie_path.exists() or cookie_path.stat().st_size == 0:
        print("Login failed — cookie file not created.")
        return False

    print(f"Login successful. Cookie saved to {cookie_path}.")
    return True
