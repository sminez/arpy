from prompt_toolkit.filters import IsDone
from prompt_toolkit.document import Document
from prompt_toolkit.history import InMemoryHistory
from prompt_toolkit.interface import CommandLineInterface
from prompt_toolkit.auto_suggest import AutoSuggestFromHistory
from prompt_toolkit.layout.processors import \
    ConditionalProcessor, HighlightMatchingBracketProcessor
from prompt_toolkit.shortcuts import \
        create_prompt_application, create_output, create_eventloop

from utils.lexparse import ArpyLexer, ArpyParser


def run(argv=""):
    '''
    The main read eval print loop for arpy.
    Uses prompt_toolkit:
        http://python-prompt-toolkit.readthedocs.io/en/stable/
    '''
    lexer = ArpyLexer()
    parser = ArpyParser()

    if argv:
        return parser.parse(lexer.tokenize(argv))
    else:
        print(' Welcome to arpy :: Calculation with Absolute Relativity\n')

    # Show matching parentheses, but only while editing.
    highlight_parens = ConditionalProcessor(
        processor=HighlightMatchingBracketProcessor(
            chars='[](){}<>'),
        filter=~IsDone())
    history = InMemoryHistory()

    while True:
        repl = create_prompt_application(
                'ξα > ',
                history=history,
                mouse_support=False,
                enable_history_search=True,
                auto_suggest=AutoSuggestFromHistory(),
                extra_input_processors=[highlight_parens])

        try:
            eventloop = create_eventloop()
            cli = CommandLineInterface(
                application=repl,
                eventloop=eventloop,
                output=create_output(true_color=True))

            user_input = cli.run(reset_current_buffer=False)
            if user_input:
                if isinstance(user_input, Document):
                    user_input = user_input.text
                lines = user_input.split('\n')
                expression = ' '.join([l.strip() for l in lines])
                parser.parse(lexer.tokenize(expression))
        except (EOFError, KeyboardInterrupt):
            # User hit Ctrl+d / Ctrl+c
            break
        finally:
            eventloop.close()

if __name__ == '__main__':
    run()
