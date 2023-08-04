from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal
from textual.message import Message
from textual.widgets import Header, Footer, Button, Label, Static

DEFAULT_INFORMATION = "Greetings from Eureka. Select any tab on the left to edit parameters for that step."

BUTTONS = {
    "Samples": "Add information about samples.",
    "Data": "Add location of data files.",
    "Basecalling": "Specify parameters for basecalling.",
    "Alignment": "Specify parameters for alignment.",
    "Consensus Calling": "Specify parameters for consensus calling.",
    "Masking": "Specify parameters for masking.",
}


class InfoButton( Button ):
    class Describe( Message ):

        def __init__( self, description: str ) -> None:
            self.description = description
            super().__init__()

    def __init__( self, label: str, description: str ) -> None:
        self.description = description
        super().__init__( label=label )

    def on_enter( self ) -> None:
        self.post_message( self.Describe( self.description ) )

    def on_leave( self ) -> None:
        self.post_message( self.Describe( DEFAULT_INFORMATION ) )


class StopwatchApp( App ):
    """A Textual app to manage stopwatches."""

    CSS_PATH = "style.css"
    TITLE = "Eureka configuration builder"
    BINDINGS = [
        Binding( "ctrl+c", "quit", "Quit", show=True, priority=True ),
        Binding( "ctrl+s", "", "Save config", show=True )
    ]

    def compose( self ) -> ComposeResult:
        """Create child widgets for the app."""
        header = Header()
        yield header
        yield Footer()

        container = Container( *(InfoButton( label, desc ) for label, desc in BUTTONS.items()), id="sidebar" )
        # container.border_title = "Tabs"
        parameters = Container( Label( "these are the parameters" ), id="parameters" )
        parameters.border_title = "Parameters"

        information = Static( DEFAULT_INFORMATION, id="information" )
        information.border_title = "Information"
        foo = Container( parameters, information )

        with Horizontal( id="mainlayout" ):
            yield container
            yield foo

    def on_info_button_describe( self, message: InfoButton.Describe ) -> None:
        print( "Called" )
        self.query_one( "#information" ).update( renderable=message.description )


if __name__ == "__main__":
    app = StopwatchApp()
    app.run()
