# REMORA Input for VS Code

A local, declarative language extension for REMORA input files. No build, npm install, or executable extension code is required.

## Install the folder

1. Extract `remora-input.zip`. The archive contains one top-level `remora-input` folder.
2. Copy that folder into your VS Code extensions directory:
   - macOS / Linux: `~/.vscode/extensions/remora-input`
   - Windows: `%USERPROFILE%\.vscode\extensions\remora-input`
   - VS Code Insiders: `~/.vscode-insiders/extensions/remora-input`
   - If you use a custom extensions directory, copy it there instead.
3. Make sure `package.json` sits directly inside that folder (avoid a second nested `remora-input` folder).
4. Restart VS Code, or run **Developer: Reload Window**.
5. Open `sample.inputs`; the status bar should show **REMORA Input**. For any other filename, click the language indicator and select **REMORA Input**.

This ZIP is a folder archive, not a VSIX; extract and copy it rather than using **Install from VSIX**. In SSH, WSL, or container windows, install the folder in the extension directory for that environment if needed.

Automatically recognized names: `inputs`, `inputs.*`, `*.inputs`, and `*.remora`. Broad `inputs.*` matching can also catch other applications' inputs files. Use the language selector to override that, or add a narrower association to your workspace's `.vscode/settings.json`:

```json
{
  "files.associations": {
    "inputs_seamount": "remora-input"
  }
}
```

## Recommended colors

TextMate scopes identify token types; your active theme decides their colors. To guarantee different colors for `remora` and `amr` and the specialized parameter groups, open **Preferences: Open User Settings (JSON)** and merge the following setting into the existing root object. If you already have `editor.tokenColorCustomizations`, append these rules to its `textMateRules` array. These colors are intended for dark themes; adjust them for light backgrounds. The rules target only this extension's scopes.

```json
{
  "editor.tokenColorCustomizations": {
    "textMateRules": [
      {
        "scope": [
          "comment.line.number-sign.remora-input"
        ],
        "settings": {
          "foreground": "#6A9955",
          "fontStyle": "italic"
        }
      },
      {
        "scope": [
          "entity.name.namespace.remora"
        ],
        "settings": {
          "foreground": "#4EC9B0",
          "fontStyle": "bold"
        }
      },
      {
        "scope": [
          "entity.name.namespace.amr.remora-input"
        ],
        "settings": {
          "foreground": "#C586C0",
          "fontStyle": "bold"
        }
      },
      {
        "scope": [
          "variable.other.parameter.remora-input"
        ],
        "settings": {
          "foreground": "#9CDCFE",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "support.type.boundary.remora-input",
          "variable.other.boundary.remora-input"
        ],
        "settings": {
          "foreground": "#DCDCAA",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "support.type.turbulence.remora-input",
          "variable.other.turbulence.remora-input"
        ],
        "settings": {
          "foreground": "#D7BA7D",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "entity.name.section.remora-input"
        ],
        "settings": {
          "foreground": "#C586C0",
          "fontStyle": "italic"
        }
      },
      {
        "scope": [
          "variable.other.nested.remora-input"
        ],
        "settings": {
          "foreground": "#9CDCFE",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "constant.numeric.remora-input"
        ],
        "settings": {
          "foreground": "#B5CEA8",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "constant.language.boolean.remora-input"
        ],
        "settings": {
          "foreground": "#569CD6",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "string.quoted.double.remora-input",
          "string.quoted.single.remora-input"
        ],
        "settings": {
          "foreground": "#CE9178",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "string.unquoted.remora-input"
        ],
        "settings": {
          "foreground": "#CE9178",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "string.unquoted.path.remora-input"
        ],
        "settings": {
          "foreground": "#D7BA7D",
          "fontStyle": ""
        }
      },
      {
        "scope": [
          "keyword.operator.assignment.remora-input",
          "punctuation.accessor.remora-input"
        ],
        "settings": {
          "foreground": "#D4D4D4",
          "fontStyle": ""
        }
      }
    ]
  }
}
```

## Highlighting behavior

- `#` starts a full-line or inline comment outside quotes. Commented-out assignments remain entirely comments.
- `remora` and `amr` have separate namespace scopes. Their dots are punctuation; parameter names have variable scopes.
- `remora.bc.*` uses boundary scopes. `remora.gls_*` uses turbulence scopes, with `gls_` separated from its suffix. Other nested names such as `remora.seamount.in_box_lo` distinguish the first section (`seamount`) from the remaining parameter name.
- Signed integers and decimals, leading/trailing decimal points, and scientific notation (`1e-6`, `-2E+3`, `1D-4`) use numeric scopes. `0` and `1` stay numeric.
- `true`, `false`, `yes`, `no`, `on`, and `off` are case-insensitive boolean tokens.
- Single- and double-quoted strings support backslash escapes. A `#` inside quotes remains string content. Unclosed quotes end at the line boundary so highlighting recovers on the next line.
- Bare keywords and values use string scopes. Unquoted paths containing `/` or `\`, and filenames with alphabetic extensions, use a separate path scope. Quoted filenames keep their quoted-string scope. Path detection is a visual heuristic, not filesystem validation.
- Generic assignments such as `geometry.coord_sys = 0` are supported too.

Assignments begin at the start of a line, allowing indentation, and use `=`. This grammar is a syntax highlighter, not a REMORA configuration validator. The sample reproduces the supplied syntax with additional synthetic cases; do not treat the whole sample as a runnable simulation. Shell-style continued value lines and multiline strings are not modeled.

Use **Developer: Inspect Editor Tokens and Scopes** to inspect the scope under the cursor. Line comment toggling uses `#`.

## Contents

- `package.json`: language registration and file associations
- `language-configuration.json`: comments, quote pairing, and word selection
- `syntaxes/remora-input.tmLanguage.json`: TextMate grammar
- `sample.inputs`: supplied-style inputs plus highlighting examples
- `README.md`: installation and color settings

Implementation follows the official [VS Code Syntax Highlight Guide](https://code.visualstudio.com/api/language-extensions/syntax-highlight-guide) and [Language Configuration Guide](https://code.visualstudio.com/api/language-extensions/language-configuration-guide).
