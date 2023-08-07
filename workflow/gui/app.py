from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal
from textual.widgets import Header, Footer, Static

from widgets import InfoButton, InfoInput
from messages import DEFAULT_INFORMATION, MEAN_MESSAGE

BUTTONS = {
    "Samples": "Add information about samples.",
    "Data": "Add location of data files.",
    "Basecalling": "Specify parameters for basecalling.",
    "Alignment": "Specify parameters for alignment.",
    "Consensus Calling": "Specify parameters for consensus calling.",
    "Masking": "Specify parameters for masking.",
}


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

        parameters = Container( InfoInput( MEAN_MESSAGE ), InfoInput( MEAN_MESSAGE ), id="parameters" )
        parameters.border_title = "Parameters"

        information = Static( DEFAULT_INFORMATION, id="information" )
        information.border_title = "Information"
        rightside = Container( parameters, information )

        with Horizontal( id="mainlayout" ):
            yield container
            yield rightside

    def on_info_button_description( self, message: InfoButton.Description ) -> None:
        self.query_one( "#information" ).update( renderable=message.description )

    def on_info_input_description( self, message: InfoInput.Description ) -> None:
        self.query_one( "#information" ).update( renderable=message.description )


if __name__ == "__main__":
    app = StopwatchApp()
    app.run()
